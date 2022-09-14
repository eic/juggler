#pragma once

#include <optional>
#include <type_traits>
#include <vector>

// Useful type traits for generically translation framework specific algorithms into
// using `algorithms'

namespace algorithms {
// Make it easy to handle the two special data types for algorithms: std::optional<T>
// and std::vector<R> where needed
template <class T> struct is_vector : std::false_type {};
template <class T, class A> struct is_vector<std::vector<T, A>> : std::true_type {};
template <class T> constexpr bool is_vector_v = is_vector<T>::value;

template <class T> struct is_optional : std::false_type {};
template <class T> struct is_optional<std::optional<T>> : std::true_type {};
template <class T> constexpr bool is_optional_v = is_optional<T>::value;

} // namespace algorithms

