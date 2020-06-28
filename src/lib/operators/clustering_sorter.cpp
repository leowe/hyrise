#include "clustering_sorter.hpp"

#include <memory>
#include <string>

#include "concurrency/transaction_context.hpp"
#include "operators/sort.hpp"
#include "operators/table_wrapper.hpp"
#include "operators/validate.hpp"
#include "storage/chunk_encoder.hpp"
#include "storage/table.hpp"
#include "utils/assert.hpp"

namespace opossum {

ClusteringSorter::ClusteringSorter(const std::shared_ptr<const AbstractOperator>& referencing_table_op, std::shared_ptr<Table> table, const std::set<ChunkID>& chunk_ids, const ColumnID sort_column_id)
    : AbstractReadWriteOperator{OperatorType::ClusteringSorter, referencing_table_op}, _table{table}, _chunk_ids{chunk_ids}, _sort_column_id{sort_column_id}, _num_locks{0}, _transaction_id{0} {
      size_t num_rows = 0;
      for (const auto chunk_id : _chunk_ids) {
        const auto& chunk = _table->get_chunk(chunk_id);
        Assert(chunk, "chunk disappeared");
        num_rows += chunk->size();
      }
    }

const std::string& ClusteringSorter::name() const {
  static const auto name = std::string{"ClusteringSorter"};
  return name;
}

std::shared_ptr<const Table> ClusteringSorter::_on_execute(std::shared_ptr<TransactionContext> context) {
  _transaction_id = context->transaction_id();

  // Sort will discard MVCC information anyway, so need to keep it
  auto sorting_table = std::make_shared<Table>(_table->column_definitions(), TableType::Data, _table->target_chunk_size(), UseMvcc::No);

  std::vector<size_t> invalid_row_counts;
  for (const auto chunk_id : _chunk_ids) {
    const auto& chunk = _table->get_chunk(chunk_id);
    Assert(chunk, "chunk must not be deleted");
    invalid_row_counts.push_back(chunk->invalid_row_count());

    Segments segments;
    for (ColumnID column_id{0}; column_id < _table->column_count(); column_id++) {
      const auto &segment = chunk->get_segment(column_id);
      segments.push_back(segment);
    }

    sorting_table->append_chunk(segments);
  }

  auto wrapper = std::make_shared<TableWrapper>(sorting_table);
  wrapper->execute();

  //auto validate = std::make_shared<Validate>(wrapper);
  //validate->execute();

  const std::vector<SortColumnDefinition> sort_column_definitions = { SortColumnDefinition(_sort_column_id, SortMode::Ascending) };
  auto sort = std::make_shared<Sort>(wrapper, sort_column_definitions, _table->target_chunk_size(), Sort::ForceMaterialization::Yes);
  sort->execute();
  _sorted_table = sort->get_output();


  // get locks for the unsorted chunks in the table
  size_t invalid_row_index = 0;
  for (const auto chunk_id : _chunk_ids) {
    const auto& chunk = _table->get_chunk(chunk_id);
    Assert(chunk, "chunk is not supposed to be deleted");

    // TODO @Jan?
    Assert(chunk->invalid_row_count() == 0, "Cannot handle invalidations in the clustered, unsorted chunks yet. This is because Sort discards MvccData (and thus invalidations). Maybe add a Validate?");

    const auto success = _lock_chunk(chunk);
    if (!success) {
      return nullptr;
    }

    if (chunk->invalid_row_count() != invalid_row_counts[invalid_row_index]) {
      // chunk was modified between sorting and locking
      _mark_as_failed();
      return nullptr;
    }
    invalid_row_index++;
  }

  // no need to get locks for the sorted chunks, as they get inserted as completely new chunks

  return nullptr;
}

void ClusteringSorter::_unlock_all() {
  // we only hold locks for the unsorted chunks
  for (const auto chunk_id : _chunk_ids) {
    const auto& chunk = _table->get_chunk(chunk_id);
    Assert(chunk, "chunk is not supposed to be deleted");
    _unlock_chunk(chunk);
  }

  Assert(_num_locks == 0, "there should be no more locks, but got " + std::to_string(_num_locks));
}

bool ClusteringSorter::_lock_chunk(const std::shared_ptr<Chunk> chunk) {
  const auto mvcc_data = chunk->mvcc_data();

  for (ChunkOffset offset{0}; offset < chunk->size(); offset++) {
    const auto expected = 0u;
    auto success = mvcc_data->compare_exchange_tid(offset, expected, _transaction_id);
    if (!success) {
      _mark_as_failed();
      return false;
    } else {
      _num_locks++;
    }
  }

  return true;
}

void ClusteringSorter::_unlock_chunk(const std::shared_ptr<Chunk> chunk) {
  const auto mvcc_data = chunk->mvcc_data();

  for (ChunkOffset offset{0}; offset < chunk->size(); offset++) {
    if (mvcc_data->get_tid(offset) == _transaction_id) {
      const auto success = mvcc_data->compare_exchange_tid(offset, _transaction_id, 0u);
      Assert(success, "Unable to unlock a row that belongs to our own transaction");
      _num_locks--;
    }
  }
}

void ClusteringSorter::_on_commit_records(const CommitID commit_id) {
  // all locks have been acquired by now

  // MVCC-delete the unsorted chunks
  for (const auto chunk_id : _chunk_ids) {
    const auto& chunk = _table->get_chunk(chunk_id);
    // TODO can this happen?
    Assert(chunk, "chunk disappeared");
    Assert(chunk->invalid_row_count() == 0, "chunk should not have invalid rows");

    const auto& mvcc_data = chunk->mvcc_data();
    for (ChunkOffset offset{0}; offset < chunk->size(); offset++) {
      mvcc_data->set_end_cid(offset, commit_id);
    }
    chunk->increase_invalid_row_count(chunk->size());
  }

  // copy the chunks from the sorted table over and update MVCC accordingly
  for (ChunkID chunk_id{0}; chunk_id < _sorted_table->chunk_count(); chunk_id++) {
    const auto& chunk = _sorted_table->get_chunk(chunk_id);
    Assert(chunk, "_sorted_table is not supposed to have removed chunks");

    Segments segments;
    for (ColumnID col_id{0}; col_id < chunk->column_count(); col_id++) {
      const auto& segment = chunk->get_segment(col_id);
      Assert(segment, "segment was null");
      segments.push_back(segment);
    }
    const auto mvcc_data = std::make_shared<MvccData>(chunk->size(), commit_id);

    // transfer meta information
    std::shared_ptr<Chunk> table_chunk;
    {
      const auto append_lock = _table->acquire_append_mutex();
      _table->append_chunk(segments, mvcc_data);
      table_chunk = _table->last_chunk();
      Assert(table_chunk, "Chunk disappeared");
    }

    table_chunk->finalize();

    Assert(!chunk->sorted_by().empty(), "chunk has no sorting information");
    table_chunk->set_sorted_by(chunk->sorted_by());

    // TODO (maybe): move encoding to disjoint_clusters_algo
    ChunkEncoder::encode_chunk(table_chunk, _table->column_data_types(), EncodingType::Dictionary);
    //Assert(chunk->pruning_statistics(), "chunk has no pruning statistics");
    //table_chunk->set_pruning_statistics(*chunk->pruning_statistics());
  }

  for (const auto chunk_id : _chunk_ids) {
    const auto& chunk = _table->get_chunk(chunk_id);
    Assert(chunk, "Chunk disappeared");
    chunk->set_cleanup_commit_id(commit_id);
  }

  _unlock_all();
}

void ClusteringSorter::_on_rollback_records() {
  _unlock_all();
}

std::shared_ptr<AbstractOperator> ClusteringSorter::_on_deep_copy(
    const std::shared_ptr<AbstractOperator>& copied_input_left,
    const std::shared_ptr<AbstractOperator>& copied_input_right) const {
  //return std::make_shared<Clustering>(copied_input_left);
  return nullptr;
}

void ClusteringSorter::_on_set_parameters(const std::unordered_map<ParameterID, AllTypeVariant>& parameters) {}

}  // namespace opossum
