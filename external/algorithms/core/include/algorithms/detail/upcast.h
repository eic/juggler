#pragma once

#include <cstdint>
#include <string>
#include <type_traits>

// Safe type conversions (without narrowing) to be used for properties,
// gets arround the limitations of type-erasure without explicitly needing
// to specify template arguments in setProperty etc.

namespace algorithms::detail {
template <typename T, typename Enable = void> struct upcast_type { using type = T; };
template <typename UInt>
struct upcast_type<UInt, std::enable_if_t<std::is_integral_v<UInt> && !std::is_signed_v<UInt>>> {
  using type = uint64_t;
};
template <typename Int>
struct upcast_type<Int, std::enable_if_t<std::is_integral_v<Int> && std::is_signed_v<Int>>> {
  using type = int64_t;
};
template <typename Enum> struct upcast_type<Enum, std::enable_if_t<std::is_enum_v<Enum>>> {
  using type = int64_t;
};
template <typename Float>
struct upcast_type<Float, std::enable_if_t<std::is_floating_point_v<Float>>> {
  using type = double;
};
template <typename String>
struct upcast_type<String, std::enable_if_t<std::is_convertible_v<String, std::string>>> {
  using type = std::string;
};
template <> struct upcast_type<bool> { using type = bool; };

template <class T> using upcast_type_t = typename upcast_type<T>::type;

template <class T> upcast_type_t<T> upcast(T&& value) { return value; }

} // namespace algorithms::detail

