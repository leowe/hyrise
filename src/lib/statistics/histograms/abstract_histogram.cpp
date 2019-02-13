#include "abstract_histogram.hpp"

#include <cmath>

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>

#include "expression/evaluation/like_matcher.hpp"
#include "generic_histogram.hpp"
#include "histogram_utils.hpp"
#include "resolve_type.hpp"
#include "single_bin_histogram.hpp"
#include "generic_histogram_builder.hpp"
#include "statistics/abstract_statistics_object.hpp"
#include "statistics/empty_statistics_object.hpp"
#include "statistics/statistics_utils.hpp"
#include "storage/create_iterable_from_segment.hpp"
#include "storage/segment_iterate.hpp"

namespace opossum {

using namespace opossum::histogram;  // NOLINT

template <typename T>
AbstractHistogram<T>::AbstractHistogram()
    : AbstractStatisticsObject(data_type_from_type<T>()) {}

template <>
AbstractHistogram<std::string>::AbstractHistogram() : AbstractStatisticsObject(DataType::String) {
  _string_domain.emplace();
}

template <>
AbstractHistogram<std::string>::AbstractHistogram(const StringHistogramDomain& string_domain)
    : AbstractStatisticsObject(DataType::String) {
  _string_domain.emplace(string_domain);
}

template <typename T>
const std::optional<StringHistogramDomain>& AbstractHistogram<T>::string_domain() const {
  return _string_domain;
}

template <typename T>
std::string AbstractHistogram<T>::description(const bool include_bin_info) const {
  std::stringstream stream;

  if constexpr (std::is_floating_point_v<T>) {
    stream.precision(20);
  }

  stream << histogram_name();
  stream << " value count: " << total_count() << ";";
  stream << " distinct count: " << total_distinct_count() << ";";
  stream << " bin count: " << bin_count() << ";";

  if (include_bin_info) {
    stream << "  Bins" << std::endl;
    for (BinID bin = 0u; bin < bin_count(); bin++) {
      if constexpr (std::is_same_v<T, std::string>) {
        stream << "  ['" << bin_minimum(bin) << "' (" << _string_domain->string_to_number(bin_minimum(bin)) << ") -> '";
        stream << bin_maximum(bin) << "' ("  << _string_domain->string_to_number(bin_maximum(bin)) << ")]: ";
      } else {
        stream << "  [" << bin_minimum(bin) << " -> " << bin_maximum(bin) << "]: ";
      }
      stream << "Height: " << bin_height(bin) << "; DistinctCount: " << bin_distinct_count(bin) << std::endl;
    }
  }

  return stream.str();
}

template <typename T>
std::vector<std::pair<T, HistogramCountType>> AbstractHistogram<T>::_gather_value_distribution(
    const std::shared_ptr<const BaseSegment>& segment, std::optional<StringHistogramDomain> string_domain) {
  std::unordered_map<T, HistogramCountType> value_counts;

  if constexpr (std::is_same_v<T, std::string>) {
    if (!string_domain) {
      string_domain.emplace();
    }
  }

  segment_iterate<T>(*segment, [&](const auto& position) {
    if (!position.is_null()) {
      if constexpr (std::is_same_v<T, std::string>) {
        ++value_counts[string_domain->string_to_domain(position.value())];
      } else {
        ++value_counts[position.value()];
      }
    }
  });

  std::vector<std::pair<T, HistogramCountType>> result(value_counts.cbegin(), value_counts.cend());

  std::sort(result.begin(), result.end(), [&](const auto& lhs, const auto& rhs) {
    return lhs.first < rhs.first;
  });

  return result;
}

template <typename T>
HistogramBin<T> AbstractHistogram<T>::bin(const BinID index) const {
  return {bin_minimum(index), bin_maximum(index), bin_height(index), bin_distinct_count(index)};
}

template <typename T>
typename AbstractHistogram<T>::HistogramWidthType AbstractHistogram<T>::bin_width(const BinID index) const {
  DebugAssert(index < bin_count(), "Index is not a valid bin.");

  // The width of an integer bin [5, 5] is 1, whereas the width of a float bin [5.1, 5.2] is 0.1

  if constexpr (std::is_same_v<T, std::string>) {
    const auto repr_min = _string_domain->string_to_number(bin_minimum(index));
    const auto repr_max = _string_domain->string_to_number(bin_maximum(index));
    return repr_max - repr_min + 1u;
  } else if constexpr (std::is_floating_point_v<T>) {
    return bin_maximum(index) - bin_minimum(index);
  } else {
    return get_next_value(bin_maximum(index) - bin_minimum(index));
  }
}

template <typename T>
T AbstractHistogram<T>::get_next_value(const T& value) const {
  if constexpr (std::is_same_v<T, std::string>) {
    return StringHistogramDomain{_string_domain->supported_characters, value.size() + 1}.next_value(value);
  } else {
    return next_value(value);
  }
}

template <typename T>
float AbstractHistogram<T>::bin_ratio_less_than(const BinID bin_id, const T& value) const {
  if (value <= bin_minimum(bin_id)) {
    return 0.0f;
  }

  if (value > bin_maximum(bin_id)) {
    return 1.0f;
  }

  if constexpr (!std::is_same_v<T, std::string>) {
    return (static_cast<float>(value) - static_cast<float>(bin_minimum(bin_id))) / static_cast<float>(bin_width(bin_id));
  } else {
   /*
    * We need to convert strings to their numerical representation to calculate a share.
    * This conversion is done based on prefixes because strings of arbitrary length cannot be converted to a numerical
    * representation that satisfies the following requirements:
    *  1. For two strings s1 and s2: s1 < s2 -> repr(s1) < repr(s2)
    *  2. For two strings s1 and s2: dist(s1, s2) == repr(s2) - repr(s1)
    *  repr(s) is the numerical representation for a string s, and dist(s1, s2) returns the number of strings between
    *  s1 and s2 in the domain of strings with at most length `string_prefix_length`
    *  and the set of supported characters `supported_characters`.
    *
    * Thus, we calculate the range based only on a domain of strings with a maximum length of `string_prefix_length`
    * characters.
    * However, we make use of a trick: if the bin edges share a common prefix, we strip that common prefix and
    * take the substring starting after that prefix.
    *
    * Example:
    *  - bin: ["intelligence", "intellij"]
    *  - supported_characters: [a-z]
    *  - string_prefix_length: 4
    *  - value: intelligent
    *
    *  Traditionally, if we did not strip the common prefix, we would calculate the range based on the
    *  substring of length `string_prefix_length`, which is "inte" for both lower and upper edge of the bin.
    *  We could not make a reasonable assumption how large the share is.
    *  Instead, we strip the common prefix ("intelli") and calculate the share based on the numerical representation
    *  of the substring after the common prefix.
    *  That is, what is the share of values smaller than "gent" in the range ["gence", "j"]?
    */

    const auto bin_min = bin_minimum(bin_id);
    const auto bin_max = bin_maximum(bin_id);

    auto common_prefix_length = size_t{0};
    const auto max_common_prefix_len = std::min(bin_min.length(), bin_max.length());
    for (; common_prefix_length < max_common_prefix_len; ++common_prefix_length) {
      if (bin_min[common_prefix_length] != bin_max[common_prefix_length]) {
        break;
      }
    }

    DebugAssert(value.substr(0, common_prefix_length) == bin_min.substr(0, common_prefix_length),
                "Value does not belong to bin");

    const auto in_domain_value = _string_domain->string_to_domain(value.substr(common_prefix_length));
    const auto value_repr = _string_domain->string_to_number(in_domain_value);
    const auto min_repr = _string_domain->string_to_number(bin_min.substr(common_prefix_length));
    const auto max_repr = _string_domain->string_to_number(bin_max.substr(common_prefix_length));
    const auto bin_ratio = static_cast<float>(value_repr - min_repr) / (max_repr - min_repr + 1);

    return bin_ratio;
  }
}

template <typename T>
float AbstractHistogram<T>::bin_ratio_less_than_equals(const BinID bin_id, const T& value) const {
  if (value < bin_minimum(bin_id)) {
    return 0.0f;
  }

  if (value >= bin_maximum(bin_id)) {
    return 1.0f;
  }

  if constexpr (!std::is_same_v<T, std::string>) {
    return bin_ratio_less_than(bin_id, next_value(value));
  } else {
    const auto bin_min = bin_minimum(bin_id);
    const auto bin_max = bin_maximum(bin_id);

    auto common_prefix_length = size_t{0};
    const auto max_common_prefix_len = std::min(bin_min.length(), bin_max.length());
    for (; common_prefix_length < max_common_prefix_len; ++common_prefix_length) {
      if (bin_min[common_prefix_length] != bin_max[common_prefix_length]) {
        break;
      }
    }

    DebugAssert(value.substr(0, common_prefix_length) == bin_min.substr(0, common_prefix_length),
                "Value does not belong to bin");

    const auto in_domain_value = _string_domain->string_to_domain(value.substr(common_prefix_length));
    const auto value_repr = _string_domain->string_to_number(in_domain_value) + 1;
    const auto min_repr = _string_domain->string_to_number(bin_min.substr(common_prefix_length));
    const auto max_repr = _string_domain->string_to_number(bin_max.substr(common_prefix_length));
    const auto bin_ratio = static_cast<float>(value_repr - min_repr) / (max_repr - min_repr + 1);

    return bin_ratio;
  }
}

template <typename T>
bool AbstractHistogram<T>::_general_does_not_contain(const PredicateCondition predicate_condition,
                                                     const AllTypeVariant& variant_value,
                                                     const std::optional<AllTypeVariant>& variant_value2) const {
  if (bin_count() == 0) {
    return true;
  }

  const auto value = type_cast_variant<T>(variant_value);

  if constexpr (std::is_same_v<T, std::string>) {
    Assert(_string_domain->contains(value), "Invalid value");
  }

  switch (predicate_condition) {
    case PredicateCondition::Equals: {
      const auto bin_id = _bin_for_value(value);
      // It is possible for EqualWidthHistograms to have empty bins.
      return bin_id == INVALID_BIN_ID || bin_height(bin_id) == 0ul;
    }
    case PredicateCondition::NotEquals:
      return bin_minimum(BinID{0}) == value && bin_maximum(bin_count() - 1) == value;
    case PredicateCondition::LessThan:
      return value <= bin_minimum(BinID{0});
    case PredicateCondition::LessThanEquals:
      return value < bin_minimum(BinID{0});
    case PredicateCondition::GreaterThanEquals:
      return value > bin_maximum(bin_count() - 1);
    case PredicateCondition::GreaterThan:
      return value >= bin_maximum(bin_count() - 1);
    case PredicateCondition::Between: {
      Assert(static_cast<bool>(variant_value2), "Between operator needs two values.");

      if (_does_not_contain(PredicateCondition::GreaterThanEquals, value)) {
        return true;
      }

      const auto value2 = type_cast_variant<T>(*variant_value2);
      if (_does_not_contain(PredicateCondition::LessThanEquals, value2) || value2 < value) {
        return true;
      }

      const auto value_bin = _bin_for_value(value);
      const auto value2_bin = _bin_for_value(value2);

      // In an EqualDistinctCountHistogram, if both values fall into the same gap, we can prune the predicate.
      // We need to have at least two bins to rule out pruning if value < min and value2 > max.
      if (value_bin == INVALID_BIN_ID && value2_bin == INVALID_BIN_ID && bin_count() > 1ul &&
          _next_bin_for_value(value) == _next_bin_for_value(value2)) {
        return true;
      }

      // In an EqualWidthHistogram, if both values fall into a bin that has no elements,
      // and there are either no bins in between or none of them have any elements, we can also prune the predicate.
      if (value_bin != INVALID_BIN_ID && value2_bin != INVALID_BIN_ID && bin_height(value_bin) == 0 &&
          bin_height(value2_bin) == 0) {
        for (auto current_bin = value_bin + 1; current_bin < value2_bin; current_bin++) {
          if (bin_height(current_bin) > 0ul) {
            return false;
          }
        }
        return true;
      }

      return false;
    }
    case PredicateCondition::Like:
    case PredicateCondition::NotLike:
      Fail("Predicate (NOT) LIKE is not supported for non-string columns.");
    default:
      // Do not prune predicates we cannot (yet) handle.
      return false;
  }
}

template <typename T>
bool AbstractHistogram<T>::_does_not_contain(const PredicateCondition predicate_condition,
                                             const AllTypeVariant& variant_value,
                                             const std::optional<AllTypeVariant>& variant_value2) const {
  return _general_does_not_contain(predicate_condition, variant_value, variant_value2);
}

template <>
bool AbstractHistogram<std::string>::_does_not_contain(const PredicateCondition predicate_condition,
                                                       const AllTypeVariant& variant_value,
                                                       const std::optional<AllTypeVariant>& variant_value2) const {
  const auto value = type_cast_variant<std::string>(variant_value);

  switch (predicate_condition) {
    case PredicateCondition::Like: {
      if (!LikeMatcher::contains_wildcard(value)) {
        return _does_not_contain(PredicateCondition::Equals, value);
      }

      // If the pattern starts with a MatchAll, we can not prune it.
      if (value.front() == '%') {
        return false;
      }

      /**
       * We can prune prefix searches iff the domain of values captured by a prefix pattern is prunable.
       *
       * Example:
       * bins: [a, b], [d, e]
       * predicate: col LIKE 'c%'
       *
       * With the same argument we can also prune predicates in the form of 'c%foo',
       * where foo can be any pattern itself.
       * We only have to consider the pattern up to the first AnyChars wildcard.
       */
      const auto match_all_index = value.find('%');
      if (match_all_index != std::string::npos) {
        const auto search_prefix = value.substr(0, std::min(match_all_index, _string_domain->prefix_length));
        if (_does_not_contain(PredicateCondition::GreaterThanEquals, search_prefix)) {
          return true;
        }

        const auto search_prefix_next_value = StringHistogramDomain{_string_domain->supported_characters, search_prefix.length()}.next_value(search_prefix);

        // If the next value is the same as the prefix, it means that there is no larger value in the domain
        // of substrings. In that case we cannot prune, because otherwise the previous check would already return true.
        if (search_prefix == search_prefix_next_value) {
          return false;
        }

        if (_does_not_contain(PredicateCondition::LessThan, search_prefix_next_value)) {
          return true;
        }

        const auto search_prefix_bin = _bin_for_value(search_prefix);
        const auto search_prefix_next_value_bin = _bin_for_value(search_prefix_next_value);

        if (search_prefix_bin == INVALID_BIN_ID) {
          const auto search_prefix_next_bin = _next_bin_for_value(search_prefix);

          // In an EqualDistinctCountHistogram, if both values fall into the same gap, we can prune the predicate.
          // We need to have at least two bins to rule out pruning if search_prefix < min
          // and search_prefix_next_value > max.
          if (search_prefix_next_value_bin == INVALID_BIN_ID && bin_count() > 1ul &&
              search_prefix_next_bin == _next_bin_for_value(search_prefix_next_value)) {
            return true;
          }

          // In an EqualDistinctCountHistogram, if the search_prefix_next_value is exactly the lower bin edge of
          // the upper bound of search_prefix, we can also prune.
          // That's because search_prefix_next_value does not belong to the range covered by the pattern,
          // but is the next value after it.
          if (search_prefix_next_value_bin != INVALID_BIN_ID &&
              search_prefix_next_bin == search_prefix_next_value_bin &&
              bin_minimum(search_prefix_next_value_bin) == search_prefix_next_value) {
            return true;
          }
        }

        // In an EqualWidthHistogram, if both values fall into a bin that has no elements,
        // and there are either no bins in between or none of them have any elements, we can also prune the predicate.
        // If the count of search_prefix_next_value_bin is not 0 but search_prefix_next_value is the lower bin edge,
        // we can still prune, because search_prefix_next_value is not part of the range (same as above).
        if (search_prefix_bin != INVALID_BIN_ID && search_prefix_next_value_bin != INVALID_BIN_ID &&
            bin_height(search_prefix_bin) == 0u &&
            (bin_height(search_prefix_next_value_bin) == 0u ||
             bin_minimum(search_prefix_next_value_bin) == search_prefix_next_value)) {
          for (auto current_bin = search_prefix_bin + 1; current_bin < search_prefix_next_value_bin; current_bin++) {
            if (bin_height(current_bin) > 0u) {
              return false;
            }
          }
          return true;
        }

        return false;
      }

      return false;
    }
    case PredicateCondition::NotLike: {
      if (!LikeMatcher::contains_wildcard(value)) {
        return _does_not_contain(PredicateCondition::NotEquals, variant_value);
      }

      // If the pattern starts with a MatchAll, we can only prune it if it matches all values.
      if (value.front() == '%') {
        return value == "%";
      }

      /**
       * We can also prune prefix searches iff the domain of values captured by the histogram is less than or equal to
       * the domain of strings captured by a prefix pattern.
       *
       * Example:
       * min: car
       * max: crime
       * predicate: col NOT LIKE 'c%'
       *
       * With the same argument we can also prune predicates in the form of 'c%foo',
       * where foo can be any pattern itself.
       * We only have to consider the pattern up to the first MatchAll character.
       */
      const auto match_all_index = value.find('%');
      if (match_all_index != std::string::npos) {
        const auto search_prefix = value.substr(0, match_all_index);
        if (search_prefix == bin_minimum(BinID{0}).substr(0, search_prefix.length()) &&
            search_prefix == bin_maximum(bin_count() - 1).substr(0, search_prefix.length())) {
          return true;
        }
      }

      return false;
    }
    default: {
      auto value_in_domain = _string_domain->string_to_domain(value);

      std::optional<AllTypeVariant> value2_in_domain;
      if (variant_value2) {
        value2_in_domain = boost::get<std::string>(*variant_value2);
      }

      return _general_does_not_contain(predicate_condition, value_in_domain, value2_in_domain);
    }
  }
}

template <typename T>
CardinalityAndDistinctCountEstimate AbstractHistogram<T>::estimate_cardinality_and_distinct_count(
    const PredicateCondition predicate_condition, const AllTypeVariant& variant_value,
    const std::optional<AllTypeVariant>& variant_value2) const {
  auto value = type_cast_variant<T>(variant_value);
  if constexpr (std::is_same_v<T, std::string>) {
    value = _string_domain->string_to_domain(value);
  }

  if (_does_not_contain(predicate_condition, variant_value, variant_value2)) {
    return {Cardinality{0}, EstimateType::MatchesNone, 0.0f};
  }

  switch (predicate_condition) {
    case PredicateCondition::Equals: {
      const auto bin_id = _bin_for_value(value);
      const auto bin_distinct_count = this->bin_distinct_count(bin_id);
    
      if (bin_distinct_count == 0) {
        return {0.0f, EstimateType::MatchesNone, 0.0f};
      } else {
        const auto cardinality = bin_height(bin_id) / bin_distinct_count;

        if (bin_distinct_count == 1) {
          return {cardinality, EstimateType::MatchesExactly, 1.0f};
        } else {
          return {cardinality, EstimateType::MatchesApproximately, 1.0f};
        }
      }
    }

    case PredicateCondition::NotEquals:
      return invert_estimate(estimate_cardinality_and_distinct_count(PredicateCondition::Equals, variant_value));

    case PredicateCondition::LessThan: {
      if (value > bin_maximum(bin_count() - 1)) {
        return {static_cast<Cardinality>(total_count()), EstimateType::MatchesAll,
                static_cast<float>(total_distinct_count())};
      }

      // This should never be false because does_not_contain should have been true further up if this was the case.
      DebugAssert(value >= bin_minimum(BinID{0}), "Value smaller than min of histogram.");

      auto cardinality = Cardinality{0};
      auto distinct_count = 0.f;
      auto estimate_type = EstimateType::MatchesApproximately;
      auto bin_id = _bin_for_value(value);

      if (bin_id == INVALID_BIN_ID) {
        // The value is within the range of the histogram, but does not belong to a bin.
        // Therefore, we need to sum up the counts of all bins with a max < value.
        bin_id = _next_bin_for_value(value);
        estimate_type = EstimateType::MatchesExactly;
      } else if (value == bin_minimum(bin_id) || bin_height(bin_id) == 0u) {
        // If the value is exactly the lower bin edge or the bin is empty,
        // we do not have to add anything of that bin and know the cardinality exactly.
        estimate_type = EstimateType::MatchesExactly;
      } else {
        const auto share = bin_ratio_less_than(bin_id, value);
        cardinality += share * bin_height(bin_id);
        distinct_count += share * bin_distinct_count(bin_id);
      }

      DebugAssert(bin_id != INVALID_BIN_ID, "Should have been caught by _does_not_contain().");

      // Sum up all bins before the bin (or gap) containing the value.
      for (BinID bin = 0u; bin < bin_id; bin++) {
        cardinality += bin_height(bin);
        distinct_count += bin_distinct_count(bin);
      }

      /**
       * The cardinality is capped at total_count().
       * It is possible for a value that is smaller than or equal to the max of the EqualHeightHistogram
       * to yield a calculated cardinality higher than total_count.
       * This is due to the way EqualHeightHistograms store the count for a bin,
       * which is in a single value (count_per_bin) for all bins rather than a vector (one value for each bin).
       * Consequently, this value is the desired count for all bins.
       * In practice, _bin_count(n) >= _count_per_bin for n < bin_count() - 1,
       * because bins are filled up until the count is at least _count_per_bin.
       * The last bin typically has a count lower than _count_per_bin.
       * Therefore, if we calculate the share of the last bin based on _count_per_bin
       * we might end up with an estimate higher than total_count(), which is then capped.
       */
      return {std::min(cardinality, static_cast<Cardinality>(total_count())), estimate_type, distinct_count};
    }
    case PredicateCondition::LessThanEquals:
      return estimate_cardinality_and_distinct_count(PredicateCondition::LessThan, get_next_value(value));

    case PredicateCondition::GreaterThanEquals:
      return invert_estimate(estimate_cardinality_and_distinct_count(PredicateCondition::LessThan, variant_value));

    case PredicateCondition::GreaterThan:
      return invert_estimate(
          estimate_cardinality_and_distinct_count(PredicateCondition::LessThanEquals, variant_value));

    case PredicateCondition::Between: {
      Assert(static_cast<bool>(variant_value2), "Between operator needs two values.");
      const auto value2 = type_cast_variant<T>(*variant_value2);

      if (value2 < value) {
        return {Cardinality{0}, EstimateType::MatchesNone, 0.0f};
      }

      const auto estimate_lt_value =
          estimate_cardinality_and_distinct_count(PredicateCondition::LessThan, variant_value);
      const auto estimate_lte_value2 =
          estimate_cardinality_and_distinct_count(PredicateCondition::LessThanEquals, *variant_value2);

      if (estimate_lt_value.type == EstimateType::MatchesAll) {
        DebugAssert(estimate_lte_value2.type == EstimateType::MatchesAll, "Estimate types do not match.");
        return {static_cast<Cardinality>(total_count()), EstimateType::MatchesAll,
                static_cast<float>(total_distinct_count())};
      }

      if (estimate_lt_value.type == EstimateType::MatchesNone) {
        DebugAssert(estimate_lte_value2.type != EstimateType::MatchesNone,
                    "Should have been caught by _does_not_contain().");
        return {estimate_lte_value2.cardinality, estimate_lte_value2.type, estimate_lte_value2.distinct_count};
      }

      if (estimate_lte_value2.type == EstimateType::MatchesNone) {
        DebugAssert(estimate_lt_value.type == EstimateType::MatchesNone, "Estimate types do not match.");
        return {Cardinality{0}, EstimateType::MatchesNone, 0.0f};
      }

      /**
       * If either estimate type is approximate, the whole result is approximate.
       * All other estimate types are exact (including MatchesNone and MatchesAll),
       * so the result estimate type is exact as well.
       */
      const auto estimate_type = estimate_lt_value.type == EstimateType::MatchesApproximately ||
                                         estimate_lte_value2.type == EstimateType::MatchesApproximately
                                     ? EstimateType::MatchesApproximately
                                     : EstimateType::MatchesExactly;

      return {estimate_lte_value2.cardinality - estimate_lt_value.cardinality, estimate_type,
              estimate_lte_value2.distinct_count - estimate_lt_value.distinct_count};
    }

    case PredicateCondition::Like:
    case PredicateCondition::NotLike:
      Fail("Predicate NOT LIKE is not supported for non-string columns.");

    default:
      Fail("Predicate not supported");
  }
}

template <typename T>
CardinalityAndDistinctCountEstimate AbstractHistogram<T>::invert_estimate(
    const CardinalityAndDistinctCountEstimate& estimate) const {
  DebugAssert(estimate.cardinality <= total_count(), "Invalid estimate cannot be inverted");
  DebugAssert(estimate.distinct_count <= total_distinct_count(), "Invalid estimate cannot be inverted");

  switch (estimate.type) {
    case EstimateType::MatchesNone:
      return {static_cast<Cardinality>(total_count()), EstimateType::MatchesAll,
              static_cast<float>(total_distinct_count())};
    case EstimateType::MatchesAll:
      return {Cardinality{0}, EstimateType::MatchesNone, 0.0f};
    case EstimateType::MatchesExactly:
    case EstimateType::MatchesApproximately:
      return {total_count() - estimate.cardinality, estimate.type, total_distinct_count() - estimate.distinct_count};
    default:
      Fail("EstimateType not supported.");
  }
}

// Specialization for numbers.
template <typename T>
CardinalityEstimate AbstractHistogram<T>::estimate_cardinality(
    const PredicateCondition predicate_condition, const AllTypeVariant& variant_value,
    const std::optional<AllTypeVariant>& variant_value2) const {
  const auto estimate = estimate_cardinality_and_distinct_count(predicate_condition, variant_value, variant_value2);
  return {estimate.cardinality, estimate.type};
}

// Specialization for strings.
template <>
CardinalityEstimate AbstractHistogram<std::string>::estimate_cardinality(
    const PredicateCondition predicate_condition, const AllTypeVariant& variant_value,
    const std::optional<AllTypeVariant>& variant_value2) const {
  const auto value = type_cast_variant<std::string>(variant_value);

  if (_does_not_contain(predicate_condition, variant_value, variant_value2)) {
    return {Cardinality{0}, EstimateType::MatchesNone};
  }

  switch (predicate_condition) {
    case PredicateCondition::Like: {
      if (!LikeMatcher::contains_wildcard(value)) {
        return estimate_cardinality(PredicateCondition::Equals, variant_value);
      }

      // We don't deal with this for now because it is not worth the effort.
      // TODO(anyone): think about good way to handle SingleChar wildcard in patterns.
      const auto single_char_count = std::count(value.cbegin(), value.cend(), '_');
      if (single_char_count > 0u) {
        return {static_cast<Cardinality>(total_count()), EstimateType::MatchesApproximately};
      }

      const auto any_chars_count = std::count(value.cbegin(), value.cend(), '%');
      DebugAssert(any_chars_count > 0u,
                  "contains_wildcard() should not return true if there is neither a '%' nor a '_' in the string.");

      // Match everything.
      if (value == "%") {
        return {static_cast<Cardinality>(total_count()), EstimateType::MatchesAll};
      }

      if (value.front() != '%') {
        /**
         * We know now we have some sort of prefix search, because there is at least one AnyChars wildcard,
         * and it is not at the start of the pattern.
         *
         * We differentiate two cases:
         *  1. Simple prefix searches, e.g., 'foo%', where there is exactly one AnyChars wildcard in the pattern,
         *  and it is at the end of the pattern.
         *  2. All others, e.g., 'foo%bar' or 'foo%bar%'.
         *
         *  The way we handle these cases is we only estimate simple prefix patterns and assume uniform distribution
         *  for additional fixed characters for the second case.
         *  Note: this is obviously far from great because not only do characters not appear with equal probability,
         *  they also appear with different probability depending on characters around them.
         *  The combination 'ing' in English is far more likely than 'qzy'.
         *  One improvement would be to have a frequency table for characters and take the probability from there,
         *  but it only gets you so far. It does not help with the second property.
         *  Nevertheless, it could be helpful especially if the number of actually occurring characters in a column are
         *  small compared to the supported characters and the frequency table would be not static but built during
         *  histogram generation.
         *  TODO(anyone): look into that in more detail.
         *
         *  That is, to estimate the first case ('foo%'), we calculate
         *  estimate_cardinality(LessThan, fop) - estimate_cardinaliy(LessThan, foo).
         *  That covers all strings starting with foo.
         *
         *  In the second case we assume that all characters in _supported_characters are equally likely to appear in
         *  a string, and therefore divide the above cardinality by the number of supported characters for each
         *  additional character that is fixed in the string after the prefix.
         *
         *  Example for 'foo%bar%baz', if we only supported the 26 lowercase latin characters:
         *  (estimate_cardinality(LessThan, fop) - estimate_cardinality(LessThan, foo)) / 26^6
         *  There are six additional fixed characters in the string ('b', 'a', 'r', 'b', 'a', and 'z').
         */
        const auto search_prefix = value.substr(0, std::min(value.find('%'), _string_domain->prefix_length));
        auto additional_characters = value.length() - search_prefix.length() - any_chars_count;

        // If there are too many fixed characters for the power to be calculated without overflow, cap the exponent.
        const auto maximum_exponent =
            std::log(std::numeric_limits<uint64_t>::max()) / std::log(_string_domain->supported_characters.length());
        if (additional_characters > maximum_exponent) {
          additional_characters = static_cast<uint64_t>(maximum_exponent);
        }

        const auto search_prefix_next_value = StringHistogramDomain{_string_domain->supported_characters, search_prefix.length()}.next_value(search_prefix);

        // If the next value is the same as the prefix, it means that there is no larger value in the domain
        // of substrings. In that case all values (total_count()) are smaller than search_prefix_next_value.
        const auto count_smaller_next_value =
            search_prefix == search_prefix_next_value
                ? total_count()
                : estimate_cardinality(PredicateCondition::LessThan, search_prefix_next_value).cardinality;

        return {
            (count_smaller_next_value - estimate_cardinality(PredicateCondition::LessThan, search_prefix).cardinality) /
                ipow(_string_domain->supported_characters.length(), additional_characters),
            EstimateType::MatchesApproximately};
      }

      /**
       * If we do not have a prefix search, but a suffix or contains search, the prefix histograms do not help us.
       * We simply assume uniform distribution for all supported characters and divide the total number of rows
       * by the number of supported characters for each additional character that is fixed (see comment above).
       *
       * Example for '%foo%b%a%', if we only supported the 26 lowercase latin characters:
       * total_count() / 26^5
       * There are five fixed characters in the string ('f', 'o', 'o', 'b', and 'a').
       */
      const auto fixed_characters = value.length() - any_chars_count;
      return {static_cast<Cardinality>(total_count()) / ipow(_string_domain->supported_characters.length(), fixed_characters),
              EstimateType::MatchesApproximately};
    }
    case PredicateCondition::NotLike: {
      if (!LikeMatcher::contains_wildcard(value)) {
        return estimate_cardinality(PredicateCondition::NotEquals, variant_value);
      }

      // We don't deal with this for now because it is not worth the effort.
      // TODO(anyone): think about good way to handle SingleChar wildcard in patterns.
      const auto single_char_count = std::count(value.cbegin(), value.cend(), '_');
      if (single_char_count > 0u) {
        return {static_cast<Cardinality>(total_count()), EstimateType::MatchesApproximately};
      }

      return invert_estimate(estimate_cardinality(PredicateCondition::Like, variant_value));
    }
    default:
      const auto estimate = estimate_cardinality_and_distinct_count(predicate_condition, variant_value, variant_value2);
      return {estimate.cardinality, estimate.type};
  }
}

template <typename T>
CardinalityEstimate AbstractHistogram<T>::invert_estimate(const CardinalityEstimate& estimate) const {
  switch (estimate.type) {
    case EstimateType::MatchesNone:
      return {static_cast<Cardinality>(total_count()), EstimateType::MatchesAll};
    case EstimateType::MatchesAll:
      return {Cardinality{0}, EstimateType::MatchesNone};
    case EstimateType::MatchesExactly:
    case EstimateType::MatchesApproximately:
      return {total_count() - estimate.cardinality, estimate.type};
    default:
      Fail("EstimateType not supported.");
  }
}

template <typename T>
std::shared_ptr<AbstractStatisticsObject> AbstractHistogram<T>::sliced(
    const PredicateCondition predicate_condition, const AllTypeVariant& variant_value,
    const std::optional<AllTypeVariant>& variant_value2) const {

  if (_does_not_contain(predicate_condition, variant_value, variant_value2)) {
    return std::make_shared<EmptyStatisticsObject>(data_type);
  }

  const auto value = type_cast_variant<T>(variant_value);

  switch (predicate_condition) {
    case PredicateCondition::Equals: {
      GenericHistogramBuilder<T> builder{1, _string_domain};
      builder.add_bin(value, value, static_cast<HistogramCountType>(
      estimate_cardinality(PredicateCondition::Equals, variant_value).cardinality), 1);
      return builder.build();
    }

    case PredicateCondition::NotEquals: {
      const auto value_bin_id = _bin_for_value(value);
      if (value_bin_id == INVALID_BIN_ID) return clone();

      const auto distinct_count = bin_distinct_count(value_bin_id);

      // Do not create empty bin, if `value` is the only value in the bin
      const auto new_bin_count = distinct_count == 1 ? bin_count() - 1 : bin_count();

      GenericHistogramBuilder<T> builder{new_bin_count, _string_domain};

      builder.add_copied_bins(*this, BinID{0}, value_bin_id);

      // Do not create empty bin.
      if (distinct_count != 1) {
        auto minimum = bin_minimum(value_bin_id);
        auto maximum = bin_maximum(value_bin_id);

        // A bin [50, 60] sliced with `!= 60` becomes [50, 59]
        // TODO(anybody) Implement bin bounds trimming for strings
        if constexpr (!std::is_same_v<std::string, T>) {
          if (minimum == value) {
            minimum = next_value(value);
          }

          if (maximum == value) {
            maximum = previous_value(value);
          }
        }

        const auto estimate = estimate_cardinality_and_distinct_count(PredicateCondition::Equals, variant_value);
        const auto new_height = bin_height(value_bin_id) - estimate.cardinality;
        const auto new_distinct_count = distinct_count - estimate.distinct_count;

        builder.add_bin(minimum, maximum, new_height, new_distinct_count);
      }

      builder.add_copied_bins(*this, value_bin_id + 1, bin_count());

      return builder.build();
    }

    case PredicateCondition::LessThanEquals:
      return sliced(PredicateCondition::LessThan, get_next_value(value));

    case PredicateCondition::LessThan: {
      auto last_bin_id = _bin_for_value(value);

      if (last_bin_id == INVALID_BIN_ID) {
        last_bin_id = _next_bin_for_value(value);

        if (last_bin_id == INVALID_BIN_ID) {
          last_bin_id = bin_count() - 1;
        } else {
          last_bin_id = last_bin_id - 1;
        }
      } 
      
      if (predicate_condition == PredicateCondition::LessThan && value == bin_minimum(last_bin_id)) {
        --last_bin_id;
      }

      auto last_bin_maximum = T{};
      // previous_value(value) is not available for strings, but we do not expect it to make a big difference.
      // TODO(anybody) Correctly implement bin bounds trimming for strings
      if constexpr (!std::is_same_v<T, std::string>) {
        last_bin_maximum = std::min(bin_maximum(last_bin_id), previous_value(value));
      } else {
        last_bin_maximum = std::min(bin_maximum(last_bin_id), value);
      }

      GenericHistogramBuilder<T> builder{last_bin_id + 1, _string_domain};
      builder.add_copied_bins(*this, BinID{0}, last_bin_id);
      builder.add_sliced_bin(*this, last_bin_id, bin_minimum(last_bin_id), last_bin_maximum);

      return builder.build();
    }

    case PredicateCondition::GreaterThan:
      return sliced(PredicateCondition::GreaterThanEquals, get_next_value(value));

    case PredicateCondition::GreaterThanEquals: {
      auto first_new_bin_id = _bin_for_value(value);

      if (first_new_bin_id == INVALID_BIN_ID) {
        first_new_bin_id = _next_bin_for_value(value);
      }

      DebugAssert(first_new_bin_id < bin_count(), "This should have been caught by _does_not_contain().");

      GenericHistogramBuilder<T> builder{bin_count() - first_new_bin_id, _string_domain};

      builder.add_sliced_bin(*this, first_new_bin_id, std::max(value, bin_minimum(first_new_bin_id)), bin_maximum(first_new_bin_id));
      builder.add_copied_bins(*this, first_new_bin_id + 1, bin_count());

      return builder.build();
    }

    case PredicateCondition::Between:
      Assert(variant_value2, "BETWEEN needs a second value.");
      return sliced(PredicateCondition::GreaterThanEquals, variant_value)
          ->sliced(PredicateCondition::LessThanEquals, *variant_value2);

    case PredicateCondition::Like:
    case PredicateCondition::NotLike:
      // TODO(anybody) Slicing for (NOT) LIKE not supported, yet
      return clone();

    case PredicateCondition::In:
    case PredicateCondition::NotIn:
    case PredicateCondition::IsNull:
    case PredicateCondition::IsNotNull:
      Fail("PredicateCondition not supported by Histograms");
  }
}

template <typename T>
std::shared_ptr<AbstractStatisticsObject> AbstractHistogram<T>::scaled(
    const Selectivity selectivity) const {
  GenericHistogramBuilder<T> builder(bin_count(), _string_domain);

  // Scale the number of values in the bin with the given selectivity.
  for (auto bin_id = BinID{0}; bin_id < bin_count(); bin_id++) {
    builder.add_bin(bin_minimum(bin_id), bin_maximum(bin_id), bin_height(bin_id) * selectivity, scale_distinct_count(selectivity, bin_height(bin_id), bin_distinct_count(bin_id)));
  }

  return builder.build();
}

template <typename T>
std::shared_ptr<AbstractHistogram<T>> AbstractHistogram<T>::split_at_bin_bounds(
    const std::vector<std::pair<T, T>>& additional_bin_edges) const {
  /**
   * Create vector with pairs for each split.
   * For a lower bin edge e, the new histogram will have to be split at previous_value(e) and e.
   * For an upper bin edge e, the new histogram will have to be split at e and next_value(e).
   *
   * We can safely ignore duplicate splits, so we use a set first and then copy the splits to a vector.
   * Also, we will create splits in gaps for now, but we will not create any actual bins for them later on.
   */
  const auto current_bin_count = bin_count();

  std::unordered_set<std::pair<T, T>, boost::hash<std::pair<T, T>>> split_set;

  for (auto bin_id = BinID{0}; bin_id < current_bin_count; bin_id++) {
    const auto bin_min = bin_minimum(bin_id);
    const auto bin_max = bin_maximum(bin_id);
    if constexpr (std::is_arithmetic_v<T>) {
      split_set.insert(std::make_pair(previous_value(bin_min), bin_min));
      split_set.insert(std::make_pair(bin_max, next_value(bin_max)));
    } else {
      // TODO(tim): turn into compile-time error
      Fail("Not supported for strings.");
    }
  }

  for (const auto& edge_pair : additional_bin_edges) {
    if constexpr (std::is_arithmetic_v<T>) {
      split_set.insert(std::make_pair(previous_value(edge_pair.first), edge_pair.first));
      split_set.insert(std::make_pair(edge_pair.second, next_value(edge_pair.second)));
    } else {
      // TODO(tim): turn into compile-time error
      Fail("Not supported for strings.");
    }
  }

  /**
   * Split pairs into single values.
   * These are the new bin edges, and we have to sort them.
   */
  std::vector<T> all_edges(split_set.size() * 2);
  size_t edge_idx = 0;
  for (const auto& split_pair : split_set) {
    all_edges[edge_idx] = split_pair.first;
    all_edges[edge_idx + 1] = split_pair.second;
    edge_idx += 2;
  }

  std::sort(all_edges.begin(), all_edges.end());

  /**
   * Remove the first and last value.
   * These are the values introduced in the first step for the first and last split.
   * The lower edge of bin 0 basically added the upper edge of bin -1,
   * and the upper edge of the last bin basically added the lower edge of the bin after the last bin.
   * Both values are obviously not needed.
   */
  all_edges.erase(all_edges.begin());
  all_edges.pop_back();


  std::vector<T> bin_minima;
  std::vector<T> bin_maxima;
  std::vector<HistogramCountType> bin_heights;
  std::vector<HistogramCountType> bin_distinct_counts;

  // We do not resize the vectors because we might not need all the slots because bins can be empty.
  const auto new_bin_count = all_edges.size() / 2;

  GenericHistogramBuilder<T> builder{new_bin_count, _string_domain};

  /**
   * Create new bins.
   * Bin edges are defined by consecutive values in pairs of two.
   */
  for (auto bin_id = BinID{0}; bin_id < new_bin_count; bin_id++) {
    const auto bin_min = all_edges[bin_id * 2];
    const auto bin_max = all_edges[bin_id * 2 + 1];

    const auto estimate = estimate_cardinality_and_distinct_count(PredicateCondition::Between, bin_min, bin_max);
    if (estimate.type == EstimateType::MatchesNone) {
      continue;
    }

    builder.add_bin(bin_min, bin_max, static_cast<HistogramCountType>(estimate.cardinality), static_cast<HistogramCountType>(estimate.distinct_count));
  }

  return builder.build();
}

template <typename T>
std::vector<std::pair<T, T>> AbstractHistogram<T>::bin_bounds() const {
  std::vector<std::pair<T, T>> bin_edges(bin_count());

  for (auto bin_id = BinID{0}; bin_id < bin_edges.size(); bin_id++) {
    bin_edges[bin_id] = std::make_pair(bin_minimum(bin_id), bin_maximum(bin_id));
  }

  return bin_edges;
}

template <typename T>
void AbstractHistogram<T>::_assert_bin_validity() {
  for (BinID bin_id{0}; bin_id < bin_count(); ++bin_id) {
    Assert(bin_minimum(bin_id) <= bin_maximum(bin_id), "Cannot have overlapping bins.");

    if (bin_id < bin_count() - 1) {
      Assert(bin_maximum(bin_id) < bin_minimum(bin_id + 1), "Bins must be sorted and cannot overlap.");
    }

    if constexpr (std::is_same_v<T, std::string>) {
      Assert(_string_domain->contains(bin_minimum(bin_id)), "Invalid string bin minimum");
      Assert(_string_domain->contains(bin_maximum(bin_id)), "Invalid string bin maximum");
    }
  }
}

EXPLICITLY_INSTANTIATE_DATA_TYPES(AbstractHistogram);

}  // namespace opossum
