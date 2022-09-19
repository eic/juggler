#pragma once

#include <array>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <algorithms/logger.h>
#include <algorithms/property.h>
#include <algorithms/type_traits.h>

namespace algorithms {

// T should either be the desired input type, a std::vector<> of the desired input type,
// or a std::optional<> of the desired input type
template <class... T> struct Input : std::tuple<input_type_t<T>...> {
  using value_type                    = std::tuple<input_type_t<T>...>;
  using data_type                     = std::tuple<T...>;
  constexpr static const size_t kSize = sizeof...(T);
};
template <class... T> struct Output : std::tuple<output_type_t<T>...> {
  using value_type                    = std::tuple<output_type_t<T>...>;
  using data_type                     = std::tuple<T...>;
  constexpr static const size_t kSize = sizeof...(T);
};

template <class Data> using DataNames = std::array<const std::string, Data::kSize>;

// TODO: C++20 Concepts version for better error handling
template <class InputType, class OutputType>
class Algorithm : public PropertyMixin, public LoggerMixin {
public:
  using InputNames  = DataNames<InputType>;
  using OutputNames = DataNames<OutputType>;

  Algorithm(std::string_view name, const InputNames& input_names, const OutputNames& output_names)
      : LoggerMixin(name), m_input_names{input_names}, m_output_names{output_names} {}

  void init();
  void process(const InputType& input, const OutputType& output);

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

