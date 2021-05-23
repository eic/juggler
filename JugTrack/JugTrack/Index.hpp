#ifndef JugTrack_Index_H
#define JugTrack_Index_H

#include <boost/container/flat_map.hpp>
#include <cstdint>

namespace Jug {

  using Index = uint32_t;

  /** Store elements that are identified by an index, e.g. in another container.
   *
   *  Each index can have zero or more associated elements. A typical case could
   *  be to store all generating particles for a hit where the hit is identified
   *  by its index in the hit container.
   */
  template <typename value_t>
  using IndexMultimap = boost::container::flat_multimap<Index, value_t>;

  /** Invert the multimap.
   *  i.e. from a -> {b...} to b -> {a...}.
   *
   *  @note This assumes that the value in the initial multimap is itself a
   *    sortable index-like object, as would be the case when mapping e.g.
   *    hit ids to particle ids/ barcodes.
   */
  template <typename value_t>
  inline boost::container::flat_multimap<value_t, Index>
  invertIndexMultimap(const IndexMultimap<value_t>& multimap)
  {
    using InverseMultimap = boost::container::flat_multimap<value_t, Index>;

    // switch key-value without enforcing the new ordering (linear copy)
    typename InverseMultimap::sequence_type unordered;
    unordered.reserve(multimap.size());
    for (auto&& [index, value] : multimap) {
      // value is now the key and the index is now the value
      unordered.emplace_back(value, index);
    }

    // adopting the unordered sequence will reestablish the correct order
    InverseMultimap inverse;
    inverse.adopt_sequence(std::move(unordered));
    return inverse;
  }

} // namespace Jug

#endif
