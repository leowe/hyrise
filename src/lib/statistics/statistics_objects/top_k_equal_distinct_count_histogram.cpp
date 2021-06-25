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
TopKEqualDistinctCountHistogram<T>::TopKEqualDistinctCountHistogram(std::shared_ptr<EqualDistinctCountHistogram<T>> histogram, 
  std::vector<T>&& top_k_names, std::vector<HistogramCountType>&& top_k_counts, const HistogramDomain<T>& domain)
    : AbstractHistogram<T>(domain) {

  _histogram = histogram;
  _top_k_names = top_k_names;
  _top_k_counts = top_k_counts;

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

  // Get the first top k values and save them into vectors
  const auto K = 10;

  std::vector<T> top_k_names(K);
  std::vector<HistogramCountType> top_k_counts(K);

  // Sort values by occurence count
  auto sorted_count_values = value_distribution;
  std::sort(sorted_count_values.begin(), sorted_count_values.end(),
            [&](const auto& l, const auto& r) { return l.second > r.second; });

  if (!sorted_count_values.empty()) {
    for(auto i = 0u; i < K; i++) {
      top_k_names[i] = sorted_count_values[i].first;
      top_k_counts[i] = sorted_count_values[i].second;

    }
  }

  // Remove Top K values from value distribution
  for (auto i = 0u; i < K; i++) {
    value_distribution.erase(std::remove(value_distribution.begin(), value_distribution.end(), std::make_pair(top_k_names[i], top_k_counts[i])), value_distribution.end());
  }

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

  auto histogram = std::make_shared<EqualDistinctCountHistogram<T>>(std::move(bin_minima), std::move(bin_maxima), 
    std::move(bin_heights), distinct_count_per_bin, bin_count_with_extra_value, domain);

  return std::make_shared<TopKEqualDistinctCountHistogram<T>>(histogram, std::move(top_k_names),std::move(top_k_counts));
}

template <typename T>
std::shared_ptr<AbstractStatisticsObject> TopKEqualDistinctCountHistogram<T>::sliced(
    const PredicateCondition predicate_condition, const AllTypeVariant& variant_value,
    const std::optional<AllTypeVariant>& variant_value2) const {
  if (AbstractHistogram<T>::does_not_contain(predicate_condition, variant_value, variant_value2)) {
    return nullptr;
  }
  const auto value = lossy_variant_cast<T>(variant_value);
  DebugAssert(value, "sliced() cannot be called with NULL");

  switch (predicate_condition) {
    case PredicateCondition::Equals: {
      GenericHistogramBuilder<T> builder{1, AbstractHistogram<T>::domain()};
      bool is_top_k_value = false; 
      for (auto i = 0u; i < _top_k_names.size(); i++) {
        if(_top_k_names[i] == value) {
          is_top_k_value = true;
          builder.add_bin(*value, *value,_top_k_counts[i], 1);
          break;
        }
      }
      if(!is_top_k_value) {
        builder.add_bin(*value, *value,
                        static_cast<HistogramCountType>(AbstractHistogram<T>::estimate_cardinality(PredicateCondition::Equals, variant_value)),
                        1);
      }
      return builder.build();
    }
    case PredicateCondition::NotEquals: {
      // Check if value in top-k values
      auto top_k_index = _top_k_names.size();
      
      for (auto i = 0u; i < _top_k_names.size(); i++) {
        if(_top_k_names[i] == value) {
          top_k_index = i;
          break;
        }
      }

      if (top_k_index != _top_k_names.size()) {
        const auto new_bin_count = bin_count() + _top_k_names.size() - 1;
        GenericHistogramBuilder<T> builder{new_bin_count, AbstractHistogram<T>::domain()};
        for (auto i = 0u, top_k_value_count = _top_k_names.size(); i < top_k_value_count; i++) {
          if (i == top_k_index) continue;
          builder.add_bin(_top_k_names[i], _top_k_names[i], _top_k_counts[i], 1);
        }
        builder.add_copied_bins(*this, BinID{0}, bin_count());
        return builder.build();
      } else {
        // Histogram or nothing

        const auto value_bin_id = _bin_for_value(*value);
        if (value_bin_id == INVALID_BIN_ID) return clone();

        auto minimum = bin_minimum(value_bin_id);
        auto maximum = bin_maximum(value_bin_id);
        const auto distinct_count = bin_distinct_count(value_bin_id);

        // Do not create empty bin if `value` is the only value in the bin
        auto new_bin_count = minimum == maximum ? bin_count() - 1 : bin_count();
        new_bin_count += _top_k_names.size();

        GenericHistogramBuilder<T> builder{new_bin_count, AbstractHistogram<T>::domain()};
        builder.add_copied_bins(*this, BinID{0}, value_bin_id);

        // Do not create empty bin if `value` is the only value in the bin
        if (minimum != maximum) {
          // A bin [50, 60] sliced with `!= 60` becomes [50, 59]
          // TODO(anybody) Implement bin bounds trimming for strings
          // NOLINTNEXTLINE clang-tidy is crazy and sees a "potentially unintended semicolon" here...
          if constexpr (!std::is_same_v<pmr_string, T>) {
            if (minimum == *value) {
              minimum = AbstractHistogram<T>::domain().next_value_clamped(*value);
            }

            if (maximum == *value) {
              maximum = AbstractHistogram<T>::domain().previous_value_clamped(*value);
            }
          }

          const auto estimate = AbstractHistogram<T>::estimate_cardinality_and_distinct_count(PredicateCondition::Equals, variant_value);
          const auto new_height = bin_height(value_bin_id) - estimate.first;
          const auto new_distinct_count = distinct_count - estimate.second;

          builder.add_bin(minimum, maximum, new_height, new_distinct_count);
        }
        builder.add_copied_bins(*this, value_bin_id + 1, bin_count());
        for (auto i = 0u, top_k_value_count = _top_k_names.size(); i < top_k_value_count; i++) {
          builder.add_bin(_top_k_names[i], _top_k_names[i], _top_k_counts[i], 1);
        }
        return builder.build();
      }
    }
    case PredicateCondition::LessThanEquals:
    case PredicateCondition::LessThan:
    case PredicateCondition::GreaterThan:
    case PredicateCondition::GreaterThanEquals:
    case PredicateCondition::BetweenInclusive:
    case PredicateCondition::BetweenLowerExclusive:
    case PredicateCondition::BetweenUpperExclusive:
    case PredicateCondition::BetweenExclusive:
    case PredicateCondition::Like:
    case PredicateCondition::NotLike:
    case PredicateCondition::In:
    case PredicateCondition::NotIn:
    case PredicateCondition::IsNull:
    case PredicateCondition::IsNotNull:
      Fail("PredicateCondition not supported by Histograms");
  }

  Fail("Invalid enum value");
}


template <typename T>
std::string TopKEqualDistinctCountHistogram<T>::name() const {
  return "TopKEqualDistinctCount";
}

template <typename T>
std::shared_ptr<AbstractHistogram<T>> TopKEqualDistinctCountHistogram<T>::clone() const {
  // The new histogram needs a copy of the data
  auto histogram_copy = _histogram->clone_as_equal_distinct_count_histogram();
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
