#pragma once

#include <cstdint>
#include <functional>
#include <random>

#include <algorithms/logger.h>
#include <algorithms/resource.h>
#include <algorithms/service.h>

namespace algorithms {

// Minimalist Random service that creates orchestrates coherent algorithm-level Generator
// creation. The Generator resources use the event Context and unique algorithm ID to
// ensure a reproducible random sequence by changing the random seed based on the Context
//
// seed = f(context, algorithmID, globalSeed)
//
// The RandomSvc exposes a global seed that can be used to move the entire seed sequence
// in a reproducible fashion.
//
// Generator itself uses the C++ 64-bit MT algorithm under the hood. The process is
// implicitly thread-safe, as each algorithm instance owns its own Generator instance.
//
class RandomSvc : public LoggedService<RandomSvc> {
public:
  void init() {
    if (m_seed.hasValue()) {
      info() << "Custom global random seed requested: " << m_seed << endmsg;
    } else {
      info("Using default random sequence.");
      setProperty("seed", 0);
    }
  }
  size_t seed() const { return m_seed; }

private:
  Property<size_t> m_seed{
      this, "seed",
      "Random seed for the internal random engine. This value does not have to be set, as a unique "
      "seed is determined from the event context, the unique algorithm ID of the calling "
      "algorithm. This seed functions as a global offset to shift the random sequence if desired. "
      "Internally, we have that seed = f(context, algorithmID, globalSeed)."};

  ALGORITHMS_DEFINE_LOGGED_SERVICE(RandomSvc)
};

// thread-safe generator front-end.
class Generator : public SvcResource<RandomSvc> {
public:
  Generator() : SvcResource<RandomSvc>{"Generator"} {}

  // Override the resource context function to change the random seed on context change
  void context(const Context& c) {
    SvcResource<RandomSvc>::context(c);
    // Subtracting the global seed is OK as worst case the unsigned integer will wrap
    // which is OK.
    m_gen.seed(c.id() - service().seed());
  }
  using SvcResource<RandomSvc>::context;

  template <class Int = int> Int uniform_int(const Int min, const Int max) const {
    std::uniform_int_distribution<Int> d{min, max};
    return d(m_gen);
  }
  template <class Float = double> Float uniform_double(const Float min, const Float max) const {
    std::uniform_real_distribution<Float> d{min, max};
    return d(m_gen);
  }
  template <class Int = int> Int poisson(const Int mean) const {
    std::poisson_distribution<Int> d{mean};
    return d(m_gen);
  }
  template <class Float = double> Float exponential(const Float lambda) const {
    std::exponential_distribution<Float> d{lambda};
    return d(m_gen);
  }
  template <class Float = double> Float gaussian(const Float mu, const Float sigma) const {
    std::normal_distribution<Float> d{mu, sigma};
    return d(m_gen);
  }

private:
  mutable std::mt19937_64 m_gen;
};

} // namespace algorithms
