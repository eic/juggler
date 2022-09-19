#pragma once

#include <gsl/gsl>
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

// Get the underlying type for each of the 3 supported cases
template <class T> struct data_type { using type = T; };
template <class T, class A> struct data_type<std::vector<T, A>> { using type = T; };
template <class T> struct data_type<std::optional<T>> { using type = T; };
template <class T> using data_type_t = typename data_type<T>::type;

// Deduce inptu and output value types
template <class T> struct deduce_type {
  using input_type  = gsl::not_null<const T*>;
  using output_type = gsl::not_null<T*>;
};
template <class T, class A> struct deduce_type<std::vector<T, A>> {
  using input_type  = const std::vector<gsl::not_null<const T*>>;
  using output_type = const std::vector<gsl::not_null<T*>>;
};
template <class T> struct deduce_type<std::optional<T>> {
  using input_type  = const T*;
  using output_type = T*;
};

template <class T> using input_type_t  = typename deduce_type<T>::input_type;
template <class T> using output_type_t = typename deduce_type<T>::output_type;

} // namespace algorithms

