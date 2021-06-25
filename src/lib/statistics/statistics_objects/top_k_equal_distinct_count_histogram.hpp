#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "abstract_histogram.hpp"
#include "types.hpp"
#include "equal_distinct_count_histogram.hpp"

namespace opossum {

class Table;

/**
 * Distinct-balanced histogram.
 * The first `bin_count_with_extra_value` contain each contain `distinct_count_per_bin + 1` distinct values,
 * all other bins contain `distinct_count_per_bin` distinct values.
 * There might be gaps between bins.
 */
template <typename T>
class TopKEqualDistinctCountHistogram : public AbstractHistogram<T> {
 public:
  using AbstractHistogram<T>::AbstractHistogram;

  TopKEqualDistinctCountHistogram(std::vector<T>&& bin_minima, std::vector<T>&& bin_maxima,
                              std::vector<HistogramCountType>&& bin_heights,
                              const HistogramCountType distinct_count_per_bin, const BinID bin_count_with_extra_value,
                              const HistogramDomain<T>& domain = {},
                              std::vector<T>&& top_k_names = {}, std::vector<HistogramCountType>&& top_k_counts = {});

  /**
   * Create an TopKEqualDistinctCountHistogram for a column (spanning all Segments) of a Table
   * @param max_bin_count   Desired number of bins. Less might be created, but never more. Must not be zero.
   */
  static std::shared_ptr<TopKEqualDistinctCountHistogram<T>> from_column(const Table& table, const ColumnID column_id,
                                                                     const BinID max_bin_count,
                                                                     const HistogramDomain<T>& domain = {});

  std::string name() const override;
  std::shared_ptr<AbstractHistogram<T>> clone() const override;
  HistogramCountType total_distinct_count() const override;
  HistogramCountType total_count() const override;

  /**
   * Returns the number of bins actually present in the histogram.
   * This number can be smaller than the number of bins requested when creating a histogram.
   * The number of bins is capped at the number of distinct values in the segment.
   * Otherwise, there would be empty bins without any benefit.
   */
  BinID bin_count() const override;

  const T& bin_minimum(const BinID index) const override;
  const T& bin_maximum(const BinID index) const override;
  HistogramCountType bin_height(const BinID index) const override;
  HistogramCountType bin_distinct_count(const BinID index) const override;

  std::shared_ptr<AbstractStatisticsObject> sliced(
    const PredicateCondition predicate_condition, const AllTypeVariant& variant_value,
    const std::optional<AllTypeVariant>& variant_value2) const override;

 protected:
  BinID _bin_for_value(const T& value) const override;
  BinID _next_bin_for_value(const T& value) const override;

 private:
  /**
   * We use multiple vectors rather than a vector of structs for ease-of-use with STL library functions.
   */

  // Min values on a per-bin basis.
  std::vector<T> _bin_minima;

  // Max values on a per-bin basis.
  std::vector<T> _bin_maxima;

  // Number of values on a per-bin basis.
  std::vector<HistogramCountType> _bin_heights;

  // Number of distinct values per bin.
  HistogramCountType _distinct_count_per_bin;

  // The first bin_count_with_extra_value bins have an additional distinct value.
  BinID _bin_count_with_extra_value;

  // Aggregated counts over all bins, to avoid redundant computation
  HistogramCountType _total_count;
  HistogramCountType _total_distinct_count;

  // Pointer to the EqualDistincCountHistogram
  // std::shared_ptr<EqualDistinctCountHistogram<T>> _histogram;

  // Vectors for top k estimation
  std::vector<T> _top_k_names;
  std::vector<HistogramCountType> _top_k_counts;

  uint16_t _k = 10;

};

}  // namespace opossum
