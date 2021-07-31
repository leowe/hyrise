#include "top_k_uniform_distribution_histogram.hpp"

#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <iterator>
#include <algorithm>

#include <tsl/robin_map.h>  // NOLINT

#include "generic_histogram.hpp"
#include "resolve_type.hpp"
#include "storage/segment_iterate.hpp"

#include "expression/evaluation/like_matcher.hpp"
#include "generic_histogram.hpp"
#include "generic_histogram_builder.hpp"
#include "lossy_cast.hpp"
#include "resolve_type.hpp"
#include "statistics/statistics_objects/abstract_statistics_object.hpp"
#include "storage/create_iterable_from_segment.hpp"
#include "storage/segment_iterate.hpp"
#include "abstract_histogram.hpp"
#include "equal_distinct_count_histogram.hpp"
#include "equal_distinct_count_histogram.cpp"

#include "../benchmark/join_order_benchmark.cpp"


namespace opossum {

template <typename T>
std::shared_ptr<GenericHistogram<T>> TopKUniformDistributionHistogram<T>::from_column(
    const Table& table, const ColumnID column_id, const HistogramDomain<T>& domain) {

  auto value_distribution = value_distribution_from_column(table, column_id, domain);

  if (value_distribution.empty()) {
    return nullptr;
  }

  // If the column holds less than K distinct values use the distinct count as TOP_K instead
  const auto k = std::min(TOP_K_CLI, value_distribution.size());

  // Get the first top k values and save them into vectors
  std::vector<T> top_k_names(k);
  std::vector<HistogramCountType> top_k_counts(k);

  // Sort values by occurrence count
  auto sorted_value_counts = value_distribution;
  std::sort(sorted_value_counts.begin(), sorted_value_counts.end(),
            [&](const auto& left_value_count, const auto& right_value_count) { return left_value_count.second > right_value_count.second; });

  // Sort TOP_K values with highest occurrence count lexicographically.
  // We later use this for more performant range predicate evaluation.
  std::sort(sorted_value_counts.begin(), sorted_value_counts.begin() + k,
            [&](const auto& left_value_count, const auto& right_value_count) { return left_value_count.first < right_value_count.first; });

  for(auto top_k_index = 0u; top_k_index < k; top_k_index++) {
    top_k_names[top_k_index] = sorted_value_counts[top_k_index].first;
    top_k_counts[top_k_index] = sorted_value_counts[top_k_index].second;
  }

  // Remove TOP_K values from value distribution
  for (auto top_k_index = 0u; top_k_index < k; top_k_index++) {
    auto value_distribution_it = remove(value_distribution.begin(), value_distribution.end(), std::make_pair(top_k_names[top_k_index], top_k_counts[top_k_index]));
    value_distribution.erase(value_distribution_it, value_distribution.end());
  }

  // Each Top K value is modeled as one bin with height as its stored count.
  // Between two Top K value bins, one bin is created for potential non Top K values using an uniform distribution assumption. 
  const auto bin_count = value_distribution.size() < 1 ? BinID{top_k_names.size()} : BinID{2*top_k_names.size() + 1};

  GenericHistogramBuilder<T> builder{bin_count, domain};
  
  const auto non_top_k_count = std::accumulate(
        value_distribution.cbegin(), 
        value_distribution.cend(), 
        HistogramCountType{0}, 
        [](HistogramCountType current_count, const std::pair<T, HistogramCountType>& value_count) { return current_count + value_count.second; });

  const auto bin_distinct_count = value_distribution.size();
  const auto count_per_non_top_k_value = bin_distinct_count != 0 ? non_top_k_count / bin_distinct_count : BinID{0};

  auto current_minimum_index = 0u;
  auto current_maximum_index = value_distribution.size() - 1;

  for (auto top_k_index = 0ul, top_k_size = top_k_names.size(); top_k_index < top_k_size; top_k_index++) {
    // current top_k value
    const auto current_top_k_value = top_k_names[top_k_index];

    // find maximum value that is still smaller than current top_k value
    auto value_dist_lower_bound = std::lower_bound(value_distribution.begin(), value_distribution.end(), current_top_k_value,
      [](const auto value_count_pair, auto value) {
        return value_count_pair.first < value;
      });

    // can't build a bin if there is no value smaller than top-k
    if (!(value_dist_lower_bound == value_distribution.begin() || std::prev(value_dist_lower_bound) - value_distribution.begin() < current_minimum_index)) {
      
      // find out how many values are before current top k value
      current_maximum_index = std::prev(value_dist_lower_bound) - value_distribution.begin();
      const auto current_distinct_values = current_maximum_index - current_minimum_index + 1;
      const auto current_bin_height = current_distinct_values * count_per_non_top_k_value;

      // add bin with values before top_k
      builder.add_bin(
        value_distribution[current_minimum_index].first, 
        value_distribution[current_maximum_index].first, 
        current_bin_height, 
        current_distinct_values
      );
    }

    // add bin for topk value
    builder.add_bin(
      current_top_k_value, 
      current_top_k_value, 
      top_k_counts[top_k_index], 
      1
    );

    // advance minimum index
    current_minimum_index = value_dist_lower_bound - value_distribution.begin();
  }

  // add last bucket if non top k values are still left after the last top k value
  if (current_minimum_index <= value_distribution.size() - 1 && value_distribution.size() > 0) {
    const auto range_maximum = value_distribution.back().first;
    const auto current_distinct_values = value_distribution.size() - 1 - current_maximum_index;
    const auto current_bin_height = current_distinct_values * count_per_non_top_k_value;
    builder.add_bin(value_distribution[current_minimum_index].first, range_maximum, current_bin_height, current_distinct_values);    
  }

  return builder.build();
}

EXPLICITLY_INSTANTIATE_DATA_TYPES(TopKUniformDistributionHistogram);

}  // namespace opossum
