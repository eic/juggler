#pragma once

#include <cstdint>
#include <functional>
#include <vector>

namespace algorithms::detail {

// Auto-refreshing cached sequence from an underlying random engine allowing for multiple instances
// to be evaluated in parallel. Specs:
//   - GenFunc is required to return N random numbers between 0 and
//     std::numeric_limits<uint_fast64_t>::max()
//   - GenFunc is responsible to deal with possible simultaneous access by multiple
//     instances of CachedBitGenerator (required to be thread-safe).
//   - The owner of a CachedGenerator instance is responsible to prevent parallel
//     calls to ::operator()
//  Implements the uniform_random_bit_generator concept
class CachedBitGenerator {
public:
  using result_type = uint_fast64_t;
  using GenFunc    = std::function<std::vector<result_type>(size_t /* N */)>;
  CachedBitGenerator(const GenFunc& gen, const size_t cache_size)
      // index starts at the end of the (empty) cache to force an immediate refresh
      // on first access
      : m_gen{gen}, m_cache(cache_size), m_index{cache_size} {}

  result_type operator()() {
    if (m_index >= m_cache.size()) {
      refresh();
    }
    return m_cache[m_index++];
  }

  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

private:
  void refresh() {
    m_cache = m_gen(m_cache.size());
    m_index = 0;
  }

  GenFunc m_gen;
  std::vector<CachedBitGenerator::result_type> m_cache;
  size_t m_index;
};

} // namespace algorithms::detail
