#pragma once

#include <cstdint>
#include <functional>
#include <mutex>
#include <random>

#include <algorithms/detail/random.h>
#include <algorithms/logger.h>
#include <algorithms/service.h>

namespace algorithms {

// Random Engine callback function:
//   - Signature: std::function<std::vector<value_type>(size_t N)> --> generates a vector
//     of N numbers
//   - RandomEngineCB is required to return N random numbers between 0 and
//     std::numeric_limits<uint_fast64_t>::max()
//   - RandomEngineCB is responsible to deal with possible simultaneous access by multiple
//     Generator instances (required to be thread-safe).
using RandomEngineCB = detail::CachedBitGenerator::GenFunc;

// thread-safe generator front-end. Requires that the underlying random engine used by
// the RandomSvc is thread-safe.
class Generator {
public:
  Generator(const RandomEngineCB& gen, const size_t cache_size) : m_gen{gen, cache_size} {}

  template <class Int = int> Int uniform_int(const Int min, const Int max) const {
    std::uniform_int_distribution<Int> d{min, max};
    std::lock_guard<std::mutex> lock{m_mutex};
    return d(m_gen);
  }
  template <class Float = double> Float uniform_double(const Float min, const Float max) const {
    std::uniform_real_distribution<Float> d{min, max};
    std::lock_guard<std::mutex> lock{m_mutex};
    return d(m_gen);
  }
  template <class Int = int> Int poisson(const Int mean) const {
    std::poisson_distribution<Int> d{mean};
    std::lock_guard<std::mutex> lock{m_mutex};
    return d(m_gen);
  }
  template <class Float = double> Float exponential(const Float lambda) const {
    std::exponential_distribution<Float> d{lambda};
    std::lock_guard<std::mutex> lock{m_mutex};
    return d(m_gen);
  }
  template <class Float = double> Float gaussian(const Float mu, const Float sigma) const {
    std::normal_distribution<Float> d{mu, sigma};
    std::lock_guard<std::mutex> lock{m_mutex};
    return d(m_gen);
  }

private:
  mutable detail::CachedBitGenerator m_gen;
  mutable std::mutex m_mutex;
};

// Random service that creates multiple Generators that are linked to a single random
// engine. The Generators are safe to be used in parallel as long as the Engine itself is
// thread-safe (this is a hard requirement for MT). The Generators avoid unnecesary locking
// by running off a auto-refreshing cached random sequence.
class RandomSvc : public LoggedService<RandomSvc> {
public:
  using value_type = detail::CachedBitGenerator::result_type;

  Generator generator() { return {m_gen, m_cache_size}; }
// FIXME fix the CMake setup so these are properly found in Gaudi
#if 0 
  void init();
  void init(const RandomEngineCB& gen);
#endif
  void init() {
    if (m_seed.hasValue()) {
      info() << "Custom random seed requested: " << m_seed << endmsg;
      m_gen = createEngine(m_seed);
    }
  }
  void init(const RandomEngineCB& gen) {
    info() << "Loading external generator function." << endmsg;
    m_gen = gen;
    if (m_seed.hasValue()) {
      warning() << "Custom random seed request ignored when using external generator function"
                << endmsg;
    }
  }

#if 0
private:
  RandomEngineCB createEngine(const size_t seed = 1);
#endif
  RandomEngineCB createEngine(const size_t seed = 1) {
    return [=](const size_t size) {
      static std::mutex m;
      static std::mt19937_64 gen{seed};
      std::lock_guard<std::mutex> lock{m};
      std::vector<value_type> ret(size);
      std::generate(ret.begin(), ret.end(), gen);
      return ret;
    };
  }
  // end of FIXME

private:
  RandomEngineCB m_gen{createEngine()};
  Property<size_t> m_seed{this, "seed", "Random seed for the internal random engine"};
  Property<size_t> m_cache_size{this, "cacheSize", 1024, "Cache size for each generator instance"};
  std::mutex m_mutex;

  ALGORITHMS_DEFINE_LOGGED_SERVICE(RandomSvc)
};

} // namespace algorithms
