// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Algorithm base class, defined as a template with a tuple of Input<> and Output<>
// parameters.
//
// Known types are:
//   - Normal data type T
//   - Optional data type std::optional<T>
//   - Vector of normal data std::vector<T>
//
// For input data, this then selects:
//   - T           --> gsl::not_null<const T*> (NOT allowed to be null)
//   - optional<T> --> const T* (allowed to be null)
//   - vector<T>   --> std::vector<gsl::not_null<const T*>> (N arguments, NOT allowed to
//                                                           be null, but can be zero
//                                                           length)
//
// Same for output data, but replace `const T*` with `T*` (mutable) everywhere.
//
// The ::process() algorithm is then provided with a tuple of both the input and the
// output pointers according to this scheme.
//
// Finally, provides provides utility traits to determine if a type Input<T...> or Output<T...> are
// an Input or Output Type (is_input_v<U> and is_output_v<U>)
//
#pragma once

#include <array>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <algorithms/logger.h>
#include <algorithms/name.h>
#include <algorithms/property.h>
#include <algorithms/type_traits.h>

namespace algorithms {

// T should either be the desired input type, a std::vector<> of the desired input type,
// or a std::optional<> of the desired input type
template <class... T> struct Input : std::tuple<input_type_t<T>...> {
  constexpr static const size_t kSize = sizeof...(T);
  using value_type                    = std::tuple<input_type_t<T>...>;
  using data_type                     = std::tuple<T...>;
  using key_type                      = std::array<const std::string, kSize>;
};
template <class... T> struct Output : std::tuple<output_type_t<T>...> {
  constexpr static const size_t kSize = sizeof...(T);
  using value_type                    = std::tuple<output_type_t<T>...>;
  using data_type                     = std::tuple<T...>;
  using key_type                      = std::array<const std::string, kSize>;
};

// TODO: C++20 Concepts version for better error handling
template <class InputType, class OutputType>
class Algorithm : public PropertyMixin, public LoggerMixin, public NameMixin {
public:
  using input_type  = InputType;
  using output_type = OutputType;
  using Input       = typename input_type::value_type;
  using Output      = typename output_type::value_type;
  using InputNames  = typename input_type::key_type;
  using OutputNames = typename output_type::key_type;

  Algorithm(std::string_view name, const InputNames& input_names, const OutputNames& output_names,
            std::string_view description)
      : LoggerMixin(name)
      , NameMixin(name, description)
      , m_input_names{input_names}
      , m_output_names{output_names} {}

  void init();
  void process(const Input& input, const Output& output);

  const InputNames& inputNames() const { return m_input_names; }
  const OutputNames& outputNames() const { return m_output_names; }

private:
  const InputNames m_input_names;
  const OutputNames m_output_names;
};

namespace detail {
  template <class T> struct is_input : std::false_type {};
  template <class... T> struct is_input<Input<T...>> : std::true_type {};
  template <class T> struct is_output : std::false_type {};
  template <class... T> struct is_output<Output<T...>> : std::true_type {};
} // namespace detail
template <class T> constexpr bool is_input_v  = detail::is_input<T>::value;
template <class T> constexpr bool is_output_v = detail::is_output<T>::value;

} // namespace algorithms

