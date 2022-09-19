#pragma once

#include <gsl/gsl>

#include <algorithms/type_traits.h>

namespace algorithms {

template <class T, class Enable = void> struct deduce_type {
  using input_type  = gsl::not_null<const T*>;
  using output_type = gsl::not_null<T*>;
};
template <class T> struct deduce_type<T, std::enable_if_t<is_vector_v<T>>> {
  using input_type  = const std::vector<gsl::not_null<const T*>>;
  using output_type = const std::vector<gsl::not_null<T*>>;
};
template <class T> struct deduce_type<T, std::enable_if_t<is_optional_v<T>>> {
  using input_type  = const T*;
  using output_type = T*;
};

template <class T> using input_type_t  = typename deduce_type<T>::input_type;
template <class T> using output_type_t = typename deduce_type<T>::output_type;

} // namespace algorithms
