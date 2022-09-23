#include <algorithms/random.h>

#include <cmath>
#include <limits>
#include <random>
#include <utility>

namespace algorithms {

void RandomSvc::init() { }

double RandomSvc::uniform(double a, double b) const {
  // initialize the random uniform number generator (runif) in a range 0 to 1
  static std::uniform_real_distribution<> runif(a, b);
  return runif(rng);
}

double RandomSvc::normal(double mu, double sigma) const {
  // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform#Implementation
  constexpr double epsilon = std::numeric_limits<double>::epsilon();
  constexpr double two_pi = 2.0 * M_PI;

  static std::uniform_real_distribution<> runif(0.0, 1.0);

  double u1, u2;
  do
  {
    u1 = runif(rng);
  }
  while (u1 <= epsilon);
  u2 = runif(rng);

  //compute z0 and z1
  auto mag = sigma * sqrt(-2.0 * log(u1));
  auto z0 = mag * cos(two_pi * u2) + mu;
  //auto z1 = mag * sin(two_pi * u2) + mu;

  return z0;
}

} // namespace algorithms

