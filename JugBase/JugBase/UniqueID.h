#ifndef JUGBASE_UNIQUE_ID
#define JUGBASE_UNIQUE_ID

// Get a unique integer identifier based on a string
// Deals with possible overflow issues
//
// Both in function and mixin form (latter is useful for algorithms)

#include <cstdint>
#include <functional>
#include <limits>
#include <string>

#include <GaudiKernel/MsgStream.h>

namespace Jug {
template <class Integer> Integer uniqueID(const std::string& s) {
  std::hash<std::string> hash_alg;
  const auto fullID = hash_alg(s);
  const Integer max = std::numeric_limits<Integer>::max();
  return static_cast<Integer>(fullID & max);
}

template <class Integer = int32_t> class AlgorithmIDMixin {
public:
  AlgorithmIDMixin(const std::string& name, MsgStream& out) : m_id{uniqueID<Integer>(name)} {
    out << "Unique ID associated with '" << name << "': " << m_id << endmsg;
  }
  AlgorithmIDMixin()                        = delete;
  AlgorithmIDMixin(const AlgorithmIDMixin&) = delete;
  AlgorithmIDMixin& operator=(const AlgorithmIDMixin&) = delete;

  Integer algorithmID() const { return m_id; }

private:
  const Integer m_id;
};
} // namespace Jug

#endif
