// Get a unique integer identifier based on a string
// Deals with possible overflow issues

#pragma once

#include <functional>
#include <limits>
#include <string>

namespace Jug {
  template <class Integer> Integer uniqueID(const std::string& s) {
    std::hash<std::string> hash_alg;
    const auto fullID = hash_alg(s);
    const Integer max = std::numeric_limits<Integer>::max();
    return static_cast<Integer>(fullID & max);
  }
}
