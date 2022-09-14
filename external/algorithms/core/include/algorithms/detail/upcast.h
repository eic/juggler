#pragma once

#include <cstdint>
#include <string>
#include <type_traits>

// Safe type conversions (without narrowing) to be used for properties,
// gets arround the limitations of type-erasure without explicitly needing
// to specify template arguments in setProperty etc.

namespace algorithms::detail {
template <typename T, typename Enable = void> struct UpcastType { using Type = T; };
template <typename UInt>
struct UpcastType<UInt, std::enable_if_t<std::is_integral_v<UInt> && !std::is_signed_v<UInt>>> {
  using Type = uint64_t;
};
template <typename Int>
struct UpcastType<Int, std::enable_if_t<std::is_integral_v<Int> && std::is_signed_v<Int>>> {
  using Type = int64_t;
};
template <typename Float>
struct UpcastType<Float, std::enable_if_t<std::is_floating_point_v<Float>>> {
  using Type = double;
};
template <typename String>
struct UpcastType<String, std::enable_if_t<std::is_convertible_v<String, std::string>>> {
  using Type = std::string;
};

template <class T> using UpcastType_t = typename UpcastType<T>::Type;

template <class T> UpcastType_t<T> upcast(const T& value) { return value; }

} // namespace algorithms::detail

