#pragma once

#include <gsl/gsl>
#include <memory>
#include <random>

#include <algorithms/logger.h>
#include <algorithms/service.h>

namespace algorithms {

class RandomSvc : public LoggedService<RandomSvc> {
public:
  // Initialize the random service
  void init();

  double uniform(double a = 0.0, double b = 1.0) const;
  double normal(double mu = 0.0, double sigma = 1.0) const;

private:
  // Standard mersenne_twister_engine seeded with rd()
  mutable std::mt19937 rng{std::random_device{}()};

  ALGORITHMS_DEFINE_LOGGED_SERVICE(RandomSvc)
};

} // namespace algorithms
