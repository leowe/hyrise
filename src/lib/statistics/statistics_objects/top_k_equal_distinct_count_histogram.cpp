#include "top_k_equal_distinct_count_histogram.hpp"

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


namespace opossum {

template <typename T>
TopKEqualDistinctCountHistogram<T>::TopKEqualDistinctCountHistogram(std::shared_ptr<AbstractHistogram<T>> histogram, 
  std::vector<T>&& top_k_names, std::vector<HistogramCountType>&& top_k_counts, const HistogramDomain<T>& domain)
    : AbstractHistogram<T>(domain), _histogram(histogram), _top_k_names(std::move(top_k_names)), _top_k_counts(std::move(top_k_counts)) {

  // TODO: moved to sliced method
  if (_histogram == nullptr) {
    GenericHistogramBuilder<T> builder{0, AbstractHistogram<T>::domain()};
    _histogram = builder.build();
  }

  const auto histogram_total_count = _histogram->total_count();
  const auto top_k_total_count = std::accumulate(_top_k_counts.cbegin(), _top_k_counts.cend(), HistogramCountType{0});
  _total_count = histogram_total_count + top_k_total_count;
  _total_distinct_count = _histogram->total_distinct_count();

}

template <typename T>
std::shared_ptr<TopKEqualDistinctCountHistogram<T>> TopKEqualDistinctCountHistogram<T>::from_column(
    const Table& table, const ColumnID column_id, const BinID max_bin_count, const HistogramDomain<T>& domain) {
  Assert(max_bin_count > 0, "max_bin_count must be greater than zero ");

  auto value_distribution = value_distribution_from_column(table, column_id, domain);

  if (value_distribution.empty()) {
    return nullptr;
  }

  // If the column holds less than K distinct values use the distinct count as TOP_K instead

  const auto k = std::min(TOP_K_DEFAULT_EQ_DIST, value_distribution.size());

  // Get the first top k values and save them into vectors
  std::vector<T> top_k_names(k);
  std::vector<HistogramCountType> top_k_counts(k);

  // Sort values by occurrence count
  auto sorted_count_values = value_distribution;
  std::sort(sorted_count_values.begin(), sorted_count_values.end(),
            [&](const auto& l, const auto& r) { return l.second > r.second; });

  // Sort TOP_K values with highest occurrence count lexicographically.
  // We later use this for more performant range predicate evaluation.
  std::sort(sorted_count_values.begin(), sorted_count_values.begin() + k,
            [&](const auto& l, const auto& r) { return l.first < r.first; });

  if (!sorted_count_values.empty()) {
    for(auto i = 0u; i < k; i++) {
      top_k_names[i] = sorted_count_values[i].first;
      top_k_counts[i] = sorted_count_values[i].second;
    }
  }

  //print top k values

  // for(auto i = 0u; i < TOP_K; i++) {
  //       std::cout << "value: " << top_k_names[i] << " with count: " << top_k_counts[i] <<std::endl;
  // }

  // Remove TOP_K values from value distribution
  for (auto i = 0u; i < k; i++) {
    auto it = remove(value_distribution.begin(), value_distribution.end(), std::make_pair(top_k_names[i], top_k_counts[i]));
    value_distribution.erase(it, value_distribution.end());
  }

  // Build equal distinct count histogram
  // If there are fewer distinct values than the number of desired bins use that instead.
  const auto bin_count =
      value_distribution.size() < max_bin_count ? static_cast<BinID>(value_distribution.size()) : max_bin_count;

  // Split values evenly among bins.
  const auto distinct_count_per_bin = static_cast<size_t>(value_distribution.size() / bin_count);
  const BinID bin_count_with_extra_value = value_distribution.size() % bin_count;

  std::vector<T> bin_minima(bin_count);
  std::vector<T> bin_maxima(bin_count);
  std::vector<HistogramCountType> bin_heights(bin_count);

  // `min_value_idx` and `max_value_idx` are indices into the sorted vector `value_distribution`
  // describing which range of distinct values goes into a bin
  auto min_value_idx = BinID{0};
  for (BinID bin_idx = 0; bin_idx < bin_count; bin_idx++) {
    auto max_value_idx = min_value_idx + distinct_count_per_bin - 1;
    if (bin_idx < bin_count_with_extra_value) {
      max_value_idx++;
    }

    // We'd like to move strings, but have to copy if we need the same string for the bin_maximum
    if (std::is_same_v<T, std::string> && min_value_idx != max_value_idx) {
      bin_minima[bin_idx] = std::move(value_distribution[min_value_idx].first);
    } else {
      bin_minima[bin_idx] = value_distribution[min_value_idx].first;
    }

    bin_maxima[bin_idx] = std::move(value_distribution[max_value_idx].first);

    bin_heights[bin_idx] =
        std::accumulate(value_distribution.cbegin() + min_value_idx, value_distribution.cbegin() + max_value_idx + 1,
                        HistogramCountType{0},
                        [](HistogramCountType a, const std::pair<T, HistogramCountType>& b) { return a + b.second; });

    min_value_idx = max_value_idx + 1;
  }

  auto histogram = std::make_shared<EqualDistinctCountHistogram<T>>(
      std::move(bin_minima), std::move(bin_maxima), std::move(bin_heights),
      static_cast<HistogramCountType>(distinct_count_per_bin), bin_count_with_extra_value);

  return std::make_shared<TopKEqualDistinctCountHistogram<T>>(histogram, std::move(top_k_names), std::move(top_k_counts));
}

template <typename T>
std::shared_ptr<AbstractStatisticsObject> TopKEqualDistinctCountHistogram<T>::sliced(
    const PredicateCondition predicate_condition, const AllTypeVariant& variant_value,
    const std::optional<AllTypeVariant>& variant_value2) const {

  // if (AbstractHistogram<T>::does_not_contain(predicate_condition, variant_value, variant_value2)) {
  //   return nullptr;
  // }

  const auto value = lossy_variant_cast<T>(variant_value);
  DebugAssert(value, "sliced() cannot be called with NULL");

  switch (predicate_condition) {
    case PredicateCondition::Equals: {
      GenericHistogramBuilder<T> builder{1, AbstractHistogram<T>::domain()};

      auto value_lower_bound = std::lower_bound(_top_k_names.begin(), _top_k_names.end(), value);
      auto value_upper_bound = std::upper_bound(_top_k_names.begin(), _top_k_names.end(), value);

      if (value_lower_bound != value_upper_bound) {
        // Value in Top-K Values
        builder.add_bin(*value, *value, _top_k_counts[value_lower_bound - _top_k_names.begin()], 1);
      } else {
        builder.add_bin(*value, *value,
                static_cast<HistogramCountType>(AbstractHistogram<T>::estimate_cardinality(PredicateCondition::Equals, variant_value)),
                1);
      }
      return builder.build();
    }
    case PredicateCondition::NotEquals: {
      // Check if value in top-k values

      auto new_top_k_names = _top_k_names;
      auto new_top_k_counts = _top_k_counts;
      auto new_histogram = _histogram->clone();

      auto value_lower_bound = std::lower_bound(new_top_k_names.begin(), new_top_k_names.end(), value);
      auto value_upper_bound = std::upper_bound(_top_k_names.begin(), _top_k_names.end(), value);

      if (value_lower_bound != value_upper_bound) {
        // Value in Top-K Values
        const auto value_index = value_lower_bound - new_top_k_names.begin();
        new_top_k_names.erase(value_lower_bound);
        new_top_k_counts.erase(new_top_k_counts.begin() + value_index);

      } else {
        // Value not in Top-K Values
        new_histogram = std::static_pointer_cast<AbstractHistogram<T>>(_histogram->sliced(predicate_condition, variant_value, variant_value2));
      }

      return std::make_shared<TopKEqualDistinctCountHistogram<T>>(new_histogram, std::move(new_top_k_names), std::move(new_top_k_counts));
    }
    case PredicateCondition::LessThanEquals: {
      // Top-K Values
      // TODO: is copy necessary?
      auto new_top_k_names = _top_k_names;
      auto new_top_k_counts = _top_k_counts;

      auto upper_bound = std::upper_bound(new_top_k_names.begin(), new_top_k_names.end(), value);
      new_top_k_names.erase(upper_bound, new_top_k_names.end());
      new_top_k_counts.resize(new_top_k_names.size());

      //print new top k values
      // std::cout << "New top k values after applying PredicateCondition::LessThanEquals" << std::endl;
      // for(auto i = 0u; i < new_top_k_names.size(); i++) {
      //   std::cout << "value: " << new_top_k_names[i] << " with count: " << new_top_k_counts[i] << std::endl;
      // }

      // Histogram Values
      auto new_histogram = std::static_pointer_cast<AbstractHistogram<T>>(_histogram->sliced(predicate_condition, variant_value, variant_value2));

      return std::make_shared<TopKEqualDistinctCountHistogram<T>>(new_histogram, std::move(new_top_k_names), std::move(new_top_k_counts));
    }
    case PredicateCondition::LessThan: {
      // Top-K Values
      auto new_top_k_names = _top_k_names;
      auto new_top_k_counts = _top_k_counts;

      auto lower_bound = std::lower_bound(new_top_k_names.begin(), new_top_k_names.end(), value);
      new_top_k_names.erase(lower_bound, new_top_k_names.end());
      new_top_k_counts.resize(new_top_k_names.size());

      //print new top k values
      // std::cout << "New top k values after applying PredicateCondition::LessThan" << std::endl;
      // for(auto i = 0u; i < new_top_k_names.size(); i++) {
      //   std::cout << "value: " << new_top_k_names[i] << " with count: " << new_top_k_counts[i] << std::endl;
      // }

      // Histogram Values
      auto new_histogram = std::static_pointer_cast<AbstractHistogram<T>>(_histogram->sliced(predicate_condition, variant_value, variant_value2));

      return std::make_shared<TopKEqualDistinctCountHistogram<T>>(new_histogram, std::move(new_top_k_names), std::move(new_top_k_counts));
    }
    case PredicateCondition::GreaterThan: {
      // Top-K Values
      auto new_top_k_names = _top_k_names;
      auto new_top_k_counts = _top_k_counts;

      auto upper_bound = std::upper_bound(new_top_k_names.begin(), new_top_k_names.end(), value);

      const auto previous_top_k_names_size = new_top_k_names.size();
      new_top_k_names.erase(new_top_k_names.begin(), upper_bound);
      const auto num_deleted_top_k_values = previous_top_k_names_size - new_top_k_names.size();
      new_top_k_counts.erase(new_top_k_counts.begin(), new_top_k_counts.begin() + num_deleted_top_k_values);

      //print new top k values
      // std::cout << "New top k values after applying PredicateCondition::GreaterThan" << std::endl;
      // for(auto i = 0u; i < new_top_k_names.size(); i++) {
      //   std::cout << "value: " << new_top_k_names[i] << " with count: " << new_top_k_counts[i] << std::endl;
      // }

      // Histogram Values
      auto new_histogram = std::static_pointer_cast<AbstractHistogram<T>>(_histogram->sliced(predicate_condition, variant_value, variant_value2));

      return std::make_shared<TopKEqualDistinctCountHistogram<T>>(new_histogram, std::move(new_top_k_names), std::move(new_top_k_counts));
    }
    case PredicateCondition::GreaterThanEquals: {
      // Top-K Values
      auto new_top_k_names = _top_k_names;
      auto new_top_k_counts = _top_k_counts;

      auto lower_bound = std::lower_bound(new_top_k_names.begin(), new_top_k_names.end(), value);

      const auto previous_top_k_names_size = new_top_k_names.size();
      new_top_k_names.erase(new_top_k_names.begin(), lower_bound);
      const auto num_deleted_top_k_values = previous_top_k_names_size - new_top_k_names.size();
      new_top_k_counts.erase(new_top_k_counts.begin(), new_top_k_counts.begin() + num_deleted_top_k_values);

      //print new top k values
      // std::cout << "New top k values after applying PredicateCondition::GreaterThanEquals" << std::endl;
      // for(auto i = 0u; i < new_top_k_names.size(); i++) {
      //   std::cout << "value: " << new_top_k_names[i] << " with count: " << new_top_k_counts[i] << std::endl;
      // }

      // Histogram Values
      auto new_histogram = std::static_pointer_cast<AbstractHistogram<T>>(_histogram->sliced(predicate_condition, variant_value, variant_value2));

      return std::make_shared<TopKEqualDistinctCountHistogram<T>>(new_histogram, std::move(new_top_k_names), std::move(new_top_k_counts));
    }
    case PredicateCondition::BetweenInclusive:
      Assert(variant_value2, "BETWEEN needs a second value.");
      return std::static_pointer_cast<TopKEqualDistinctCountHistogram<T>>(sliced(PredicateCondition::GreaterThanEquals, variant_value, variant_value2))
          ->sliced(PredicateCondition::LessThanEquals, *variant_value2, variant_value2);
    case PredicateCondition::BetweenLowerExclusive:
      Assert(variant_value2, "BETWEEN needs a second value.");
      return std::static_pointer_cast<TopKEqualDistinctCountHistogram<T>>(sliced(PredicateCondition::GreaterThan, variant_value, variant_value2))
          ->sliced(PredicateCondition::LessThanEquals, *variant_value2, variant_value2);
    case PredicateCondition::BetweenUpperExclusive:
      Assert(variant_value2, "BETWEEN needs a second value.");
      return std::static_pointer_cast<TopKEqualDistinctCountHistogram<T>>(sliced(PredicateCondition::GreaterThanEquals, variant_value, variant_value2))
          ->sliced(PredicateCondition::LessThan, *variant_value2, variant_value2);
    case PredicateCondition::BetweenExclusive:
      Assert(variant_value2, "BETWEEN needs a second value.");
      return std::static_pointer_cast<TopKEqualDistinctCountHistogram<T>>(sliced(PredicateCondition::GreaterThan, variant_value, variant_value2))
          ->sliced(PredicateCondition::LessThan, *variant_value2, variant_value2);
    case PredicateCondition::Like:
    case PredicateCondition::NotLike:
      // TODO(anybody) Slicing for (NOT) LIKE not supported, yet
      return clone();
    case PredicateCondition::In:
    case PredicateCondition::NotIn:
    case PredicateCondition::IsNull:
    case PredicateCondition::IsNotNull:
      Fail("PredicateCondition not supported by TopKUnifromDistributionHistogram");
  }
  
  Fail("Invalid enum value");
}

template <typename T>
std::shared_ptr<AbstractStatisticsObject> TopKEqualDistinctCountHistogram<T>::scaled(const Selectivity selectivity) const {
  Assert(!std::isnan(selectivity), "Selectivity should not be NaN");

  // As we have no better information we assume that the occurrences of top-k values are scaled uniformly.
  auto scaled_top_k_names = _top_k_names;
  auto scaled_top_k_counts = _top_k_counts;
  std::transform(_top_k_counts.begin(), _top_k_counts.end(), scaled_top_k_counts.begin(),
                   [&selectivity](HistogramCountType count) -> HistogramCountType { return count * selectivity; });

  // Scale histogram values
  auto scaled_histogram = std::static_pointer_cast<AbstractHistogram<T>>(_histogram->scaled(selectivity));

  return std::make_shared<TopKEqualDistinctCountHistogram<T>>(scaled_histogram, std::move(scaled_top_k_names), std::move(scaled_top_k_counts));
}

template <typename T>
std::shared_ptr<AbstractStatisticsObject> TopKEqualDistinctCountHistogram<T>::pruned(
    const size_t num_values_pruned, const PredicateCondition predicate_condition, const AllTypeVariant& variant_value,
    const std::optional<AllTypeVariant>& variant_value2) const {

  // for now, we don't prune the top_k values
  auto pruned_top_k_names = _top_k_names;
  auto pruned_top_k_counts = _top_k_counts;

  auto pruned_histogram = std::static_pointer_cast<AbstractHistogram<T>>(_histogram->pruned(num_values_pruned, predicate_condition, variant_value, variant_value2));
  return std::make_shared<TopKEqualDistinctCountHistogram<T>>(pruned_histogram, std::move(pruned_top_k_names), std::move(pruned_top_k_counts));
}


template <typename T>
std::string TopKEqualDistinctCountHistogram<T>::name() const {
  return "TopKEqualDistinctCount";
}

template <typename T>
std::shared_ptr<AbstractHistogram<T>> TopKEqualDistinctCountHistogram<T>::clone() const {
  // The new histogram needs a copy of the data
  auto histogram_copy = std::static_pointer_cast<AbstractHistogram<T>>(_histogram->clone());
  auto top_k_names_copy = _top_k_names;
  auto top_k_counts_copy = _top_k_counts;
  //return histogram_copy;
  return std::make_shared<TopKEqualDistinctCountHistogram<T>>(histogram_copy, std::move(top_k_names_copy), std::move(top_k_counts_copy));
}

template <typename T>
BinID TopKEqualDistinctCountHistogram<T>::bin_count() const {
  return _histogram->bin_count();
}

template <typename T>
BinID TopKEqualDistinctCountHistogram<T>::_bin_for_value(const T& value) const {
  //TODO: look at where is that used
  for (auto i = 0u; i < _top_k_names.size(); i++) {
    if(_top_k_names[i] == value) {
      return 1;
    }
  }
  return _histogram->bin_for_value(value);
}

template <typename T>
BinID TopKEqualDistinctCountHistogram<T>::_next_bin_for_value(const T& value) const {
  //TODO: look at where is that used
  for (auto i = 0u; i < _top_k_names.size(); i++) {
    if(_top_k_names[i] == value) {
      return 1;
    }
  }
  return _histogram->next_bin_for_value(value);
}

template <typename T>
const T& TopKEqualDistinctCountHistogram<T>::bin_minimum(const BinID index) const {
  return _histogram->bin_minimum(index);
}

template <typename T>
const T& TopKEqualDistinctCountHistogram<T>::bin_maximum(const BinID index) const {
  return _histogram->bin_maximum(index);
}

template <typename T>
HistogramCountType TopKEqualDistinctCountHistogram<T>::bin_height(const BinID index) const {
  return _histogram->bin_height(index);
}

template <typename T>
HistogramCountType TopKEqualDistinctCountHistogram<T>::bin_distinct_count(const BinID index) const {
  return _histogram->bin_distinct_count(index);
}

template <typename T>
HistogramCountType TopKEqualDistinctCountHistogram<T>::total_count() const {
  return _total_count;
}

template <typename T>
HistogramCountType TopKEqualDistinctCountHistogram<T>::total_distinct_count() const {
  return _total_distinct_count;
}

EXPLICITLY_INSTANTIATE_DATA_TYPES(TopKEqualDistinctCountHistogram);

}  // namespace opossum
