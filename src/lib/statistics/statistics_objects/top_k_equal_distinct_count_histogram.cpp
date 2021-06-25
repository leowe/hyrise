#include "top_k_equal_distinct_count_histogram.hpp"

#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

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
TopKEqualDistinctCountHistogram<T>::TopKEqualDistinctCountHistogram(std::vector<T>&& bin_minima, std::vector<T>&& bin_maxima,
                                                            std::vector<HistogramCountType>&& bin_heights,
                                                            const HistogramCountType distinct_count_per_bin,
                                                            const BinID bin_count_with_extra_value,
                                                            const HistogramDomain<T>& domain,
                                                            std::vector<T>&& top_k_names, std::vector<HistogramCountType>&& top_k_counts)
    : AbstractHistogram<T>(domain),
      _bin_minima(std::move(bin_minima)),
      _bin_maxima(std::move(bin_maxima)),
      _bin_heights(std::move(bin_heights)),
      _distinct_count_per_bin(distinct_count_per_bin),
      _bin_count_with_extra_value(bin_count_with_extra_value) {
  Assert(_bin_minima.size() == _bin_maxima.size(), "Must have the same number of lower as upper bin edges.");
  Assert(_bin_minima.size() == _bin_heights.size(), "Must have the same number of edges and heights.");
  Assert(_distinct_count_per_bin > 0, "Cannot have bins with no distinct values.");
  Assert(_bin_count_with_extra_value < _bin_minima.size(), "Cannot have more bins with extra value than bins.");

  AbstractHistogram<T>::_assert_bin_validity();

  // std::vector<T> bin_minima_vector;
  // std::copy(_bin_minima.begin(), _bin_minima.end(),
  //             std::back_inserter(bin_minima_vector));

  //  std::vector<T> bin_maxima_vector;
  // std::copy(_bin_maxima.begin(), _bin_maxima.end(),
  //             std::back_inserter(bin_maxima_vector));

  //  std::vector<T> bin_heights_vector;
  // std::copy(_bin_heights.begin(), _bin_heights.end(),
  //             std::back_inserter(bin_heights_vector));

  // _histogram = std::make_shared<EqualDistinctCountHistogram<T>>(std::move(bin_minima_vector), std::move(bin_maxima_vector), 
  //   std::move(bin_heights_vector), _distinct_count_per_bin, _bin_count_with_extra_value, domain);

  _total_count = std::accumulate(_bin_heights.cbegin(), _bin_heights.cend(), HistogramCountType{0});
  _total_distinct_count = _distinct_count_per_bin * static_cast<HistogramCountType>(bin_count()) +
                          static_cast<HistogramCountType>(_bin_count_with_extra_value);

  _top_k_names = top_k_names;
  _top_k_counts = top_k_counts;
}

template <typename T>
std::shared_ptr<TopKEqualDistinctCountHistogram<T>> TopKEqualDistinctCountHistogram<T>::from_column(
    const Table& table, const ColumnID column_id, const BinID max_bin_count, const HistogramDomain<T>& domain) {
  Assert(max_bin_count > 0, "max_bin_count must be greater than zero ");

  const auto value_distribution = value_distribution_from_column(table, column_id, domain);

  const auto K = 10;

  std::vector<T> top_k_names(K);
  std::vector<HistogramCountType> top_k_counts(K);

  auto sorted_count_values = value_distribution;
  std::sort(sorted_count_values.begin(), sorted_count_values.end(),
            [&](const auto& l, const auto& r) { return l.second > r.second; });

  if (!sorted_count_values.empty()) {
    for(auto i = 0u; i < K; i++) {
      top_k_names[i] = sorted_count_values[i].first;
      top_k_counts[i] = sorted_count_values[i].second;
    }
  } else {
    std::cout << "Empty" << std::endl;
  }

  if (value_distribution.empty()) {
    return nullptr;
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

  return std::make_shared<TopKEqualDistinctCountHistogram<T>>(
      std::move(bin_minima), std::move(bin_maxima), std::move(bin_heights),
      static_cast<HistogramCountType>(distinct_count_per_bin), bin_count_with_extra_value,
      domain, std::move(top_k_names),std::move(top_k_counts));
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
    case PredicateCondition::NotEquals:
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
  auto bin_minima_copy = _bin_minima;
  auto bin_maxima_copy = _bin_maxima;
  auto bin_heights_copy = _bin_heights;

  return std::make_shared<TopKEqualDistinctCountHistogram<T>>(std::move(bin_minima_copy), std::move(bin_maxima_copy),
                                                          std::move(bin_heights_copy), _distinct_count_per_bin,
                                                          _bin_count_with_extra_value);
}

template <typename T>
BinID TopKEqualDistinctCountHistogram<T>::bin_count() const {
  return _bin_heights.size();
}

template <typename T>
BinID TopKEqualDistinctCountHistogram<T>::_bin_for_value(const T& value) const {
  const auto it = std::lower_bound(_bin_maxima.cbegin(), _bin_maxima.cend(), value);
  const auto index = static_cast<BinID>(std::distance(_bin_maxima.cbegin(), it));

  if (it == _bin_maxima.cend() || value < bin_minimum(index) || value > bin_maximum(index)) {
    return INVALID_BIN_ID;
  }

  return index;
}

template <typename T>
BinID TopKEqualDistinctCountHistogram<T>::_next_bin_for_value(const T& value) const {
  const auto it = std::upper_bound(_bin_maxima.cbegin(), _bin_maxima.cend(), value);

  if (it == _bin_maxima.cend()) {
    return INVALID_BIN_ID;
  }

  return static_cast<BinID>(std::distance(_bin_maxima.cbegin(), it));
}

template <typename T>
const T& TopKEqualDistinctCountHistogram<T>::bin_minimum(const BinID index) const {
  DebugAssert(index < _bin_minima.size(), "Index is not a valid bin.");
  return _bin_minima[index];
}

template <typename T>
const T& TopKEqualDistinctCountHistogram<T>::bin_maximum(const BinID index) const {
  DebugAssert(index < _bin_maxima.size(), "Index is not a valid bin.");
  return _bin_maxima[index];
}

template <typename T>
HistogramCountType TopKEqualDistinctCountHistogram<T>::bin_height(const BinID index) const {
  DebugAssert(index < _bin_heights.size(), "Index is not a valid bin.");
  return _bin_heights[index];
}

template <typename T>
HistogramCountType TopKEqualDistinctCountHistogram<T>::bin_distinct_count(const BinID index) const {
  DebugAssert(index < bin_count(), "Index is not a valid bin.");
  return HistogramCountType{_distinct_count_per_bin + (index < _bin_count_with_extra_value ? 1 : 0)};
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
