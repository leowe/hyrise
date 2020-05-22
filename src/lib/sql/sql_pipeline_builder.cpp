#include "sql_pipeline_builder.hpp"
#include "hyrise.hpp"
#include "utils/tracing/probes.hpp"

namespace opossum {

SQLPipelineBuilder::SQLPipelineBuilder(const std::string& sql)
    : _sql(sql), _pqp_cache(Hyrise::get().default_pqp_cache), _lqp_cache(Hyrise::get().default_lqp_cache) {}

SQLPipelineBuilder& SQLPipelineBuilder::with_mvcc(const UseMvcc use_mvcc) {
  _use_mvcc = use_mvcc;
  return *this;
}

SQLPipelineBuilder& SQLPipelineBuilder::with_optimizer(const std::shared_ptr<Optimizer>& optimizer) {
  _optimizer = optimizer;
  return *this;
}

SQLPipelineBuilder& SQLPipelineBuilder::with_transaction_context(
    const std::shared_ptr<TransactionContext>& transaction_context) {
  _transaction_context = transaction_context;
  _use_mvcc = UseMvcc::Yes;

  return *this;
}

SQLPipelineBuilder& SQLPipelineBuilder::with_pqp_cache(const std::shared_ptr<SQLPhysicalPlanCache>& pqp_cache) {
  _pqp_cache = pqp_cache;
  return *this;
}

SQLPipelineBuilder& SQLPipelineBuilder::with_lqp_cache(const std::shared_ptr<SQLLogicalPlanCache>& lqp_cache) {
  _lqp_cache = lqp_cache;
  return *this;
}

SQLPipelineBuilder& SQLPipelineBuilder::disable_mvcc() { return with_mvcc(UseMvcc::No); }

SQLPipeline SQLPipelineBuilder::create_pipeline() const {
  DTRACE_PROBE1(HYRISE, CREATE_PIPELINE, reinterpret_cast<uintptr_t>(this));
  auto optimizer = _optimizer ? _optimizer : Optimizer::create_default_optimizer();
  auto post_caching_optimizer =
      _post_caching_optimizer ? _post_caching_optimizer : Optimizer::create_post_caching_optimizer();
  auto pipeline =
      SQLPipeline(_sql, _transaction_context, _use_mvcc, optimizer, post_caching_optimizer, _pqp_cache, _lqp_cache);
  DTRACE_PROBE3(HYRISE, PIPELINE_CREATION_DONE, pipeline.get_sql_per_statement().size(), _sql.c_str(),
                reinterpret_cast<uintptr_t>(this));
  return pipeline;
}

SQLPipelineStatement SQLPipelineBuilder::create_pipeline_statement(
    std::shared_ptr<hsql::SQLParserResult> parsed_sql) const {
  auto optimizer = _optimizer ? _optimizer : Optimizer::create_default_optimizer();
  auto post_caching_optimizer =
      _post_caching_optimizer ? _post_caching_optimizer : Optimizer::create_post_caching_optimizer();

  return {_sql,      std::move(parsed_sql),  _use_mvcc,  _transaction_context,
          optimizer, post_caching_optimizer, _pqp_cache, _lqp_cache};
}

}  // namespace opossum
