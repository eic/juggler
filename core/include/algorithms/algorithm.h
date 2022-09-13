#pragma once

#include <array>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <algorithms/logger.h>
#include <algorithms/property.h>

namespace algorithms {

// T should either be the desired input type, a std::vector<> of the desired input type,
// or a std::optional<> of the desired input type
template <class... T> struct Input : std::tuple<const T&...> {
  constexpr static const size_t kSize = sizeof...(T);
  using Type                          = std::tuple<const T&...>;
  using ValueType                     = std::tuple<T...>;
};
template <class... T> struct Output : std::tuple<T&...> {
  constexpr static const size_t kSize = sizeof...(T);
  using Type                          = std::tuple<T&...>;
  using ValueType                     = std::tuple<T...>;
};

// TODO: C++20 Concepts version for better error handling
template <class InputType, class OutputType>
class Algorithm : public PropertyMixin, public LoggerMixin {
public:
  constexpr static const size_t kInputLength  = InputType::kSize;
  constexpr static const size_t kOutputLength = OutputType::kSize;

  using InputNames  = std::array<const std::string, kInputLength>;
  using OutputNames = std::array<const std::string, kOutputLength>;

  Algorithm(std::string_view name, const InputNames& input_names, const OutputNames& output_names)
      : LoggerMixin(name), m_input_names{input_names}, m_output_names{output_names} {}

  void init();
  void process(const InputType& input, OutputType& output);

  const InputNames& inputNames() const { return m_input_names; }
  const OutputNames& outputNames() const { return m_output_names; }

private:
  const InputNames m_input_names;
  const OutputNames m_output_names;
};
} // namespace algorithms

