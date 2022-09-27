#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <gsl/gsl>
#include <mutex>
#include <random>
#include <vector>

#include <algorithms/logger.h>
#include <algorithms/service.h>

namespace algorithms {

// Cached link to the underlying generator allowing for multiple instances to be evaluated
// in parallel. Specs:
//   - GenFunc is required to return N random numbers between 0 and
//     std::numeric_limits<uint_fast64_t>::max()
//   - GenFunc is responsible to deal with possible simultaneous access by multiple
//     instances of CachedBitGenerator (required to be thread-safe).
//   - The owner of a CachedGenerator instance is responsible to prevent parallel
//     calls to ::operator()
//  Implements the uniform_random_bit_generator concept
class CachedBitGenerator {
public:
  using value_type = uint_fast64_t;
  using GenFunc    = std::function<std::vector<value_type>(size_t /* N */)>;
  CachedBitGenerator(const GenFunc& gen, const size_t cache_size)
      // index starts at the end of the (empty) cache to force an immediate refresh
      // on first access
      : m_gen{gen}, m_cache(cache_size), m_index{cache_size} {}

  value_type operator()() {
    if (m_index >= m_cache.size()) {
      refresh();
    }
    return m_cache[m_index++];
  }

  constexpr value_type min() const { return 0; }
  constexpr value_type max() const { return std::numeric_limits<value_type>::max(); }

private:
  void refresh() {
    m_cache = m_gen(m_cache.size());
    m_index = 0;
  }

  GenFunc m_gen;
  std::vector<CachedBitGenerator::value_type> m_cache;
  size_t m_index;
};

// thread-safe random generator
class RandomSvc : public LoggedService<RandomSvc> {
public:
  using value_type = CachedBitGenerator::value_type;
  using GenFunc    = CachedBitGenerator::GenFunc;

  void init();
  void init(const GenFunc& gen);

public:
  class Generator {
  public:
    Generator(const GenFunc& gen, const size_t cache_size) : m_gen{gen, cache_size} {}

    template <class Int = int> Int uniform_int(const Int min, const Int max) {
      std::uniform_int_distribution<Int> d{min, max};
      std::lock_guard<std::mutex> lock{m_mutex};
      return d(m_gen);
    }
    template <class Float = double> Float uniform_double(const Float min, const Float max) {
      std::uniform_real_distribution<Float> d{min, max};
      std::lock_guard<std::mutex> lock{m_mutex};
      return d(m_gen);
    }
    template <class Int = int> Int poisson(const Int mean) {
      std::poisson_distribution<Int> d{mean};
      std::lock_guard<std::mutex> lock{m_mutex};
      return d(m_gen);
    }
    template <class Float = double> Float exponential(const Float lambda) {
      std::exponential_distribution<Float> d{lambda};
      std::lock_guard<std::mutex> lock{m_mutex};
      return d(m_gen);
    }
    template <class Float = double> Float gaussian(const Float mu, const Float sigma) {
      std::normal_distribution<Float> d{mu, sigma};
      std::lock_guard<std::mutex> lock{m_mutex};
      return d(m_gen);
    }

  private:
    CachedBitGenerator m_gen;
    std::mutex m_mutex;
  };
  Generator generator() { return {m_gen, m_cache_size}; }

private:
  GenFunc createBitGenerator(const size_t seed = 0);

  GenFunc m_gen{createBitGenerator()};
  Property<size_t> m_seed{this, "seed"};
  Property<size_t> m_cache_size{this, "cacheSize", 1024};
  std::mutex m_mutex;

  ALGORITHMS_DEFINE_LOGGED_SERVICE(RandomSvc)
};

} // namespace algorithms
