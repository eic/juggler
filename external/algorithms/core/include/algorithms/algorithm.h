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
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <tuple>
#include <typeindex>
#include <typeinfo>
#include <vector>

#include <algorithms/detail/demangle.h>
#include <algorithms/logger.h>
#include <algorithms/name.h>
#include <algorithms/property.h>
#include <algorithms/service.h>
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

class AlgorithmBase : public PropertyMixin, public LoggerMixin, public NameMixin {
public:
  AlgorithmBase(std::string_view name, std::string_view description)
      : LoggerMixin(name), NameMixin(name, description) {}
};

// TODO: C++20 Concepts version for better error handling
template <class InputType, class OutputType> class Algorithm : public AlgorithmBase {
public:
  using input_type     = InputType;
  using output_type    = OutputType;
  using algorithm_type = Algorithm;
  using Input          = typename input_type::value_type;
  using Output         = typename output_type::value_type;
  using InputNames     = typename input_type::key_type;
  using OutputNames    = typename output_type::key_type;

  Algorithm(std::string_view name, const InputNames& input_names, const OutputNames& output_names,
            std::string_view description)
      : AlgorithmBase(name, description)
      , m_input_names{input_names}
      , m_output_names{output_names} {}

  virtual ~Algorithm() {}
  virtual void init() {}
  virtual void process(const Input&, const Output&) const {}

  const InputNames& inputNames() const { return m_input_names; }
  const OutputNames& outputNames() const { return m_output_names; }

private:
  const InputNames m_input_names;
  const OutputNames m_output_names;
};

// Algorithm service that stores factories for different algorithms indexed by their
// underlying algorithm_type. This allows a framework to only need to know about every
// algorithm_type (algorithm signature) once, which is sufficient to constrain all
// compile-time code needed.
class AlgorithmSvc : public LoggedService<AlgorithmSvc> {
public:
  // Add a new factory for an algorithm with a specified type. The full demangled type name
  // will be used as identifier for the factory
  template <class Algo> void add() {
    static std::mutex m;
    static std::lock_guard<std::mutex> lock{m};
    auto type            = Algo::algorithm_type;
    std::type_index name = detail::demangledName<Algo>();
    if (!available<type>()) {
      m_factories[type] = {};
    }
    if (available<type>(name)) {
      return; // do nothing, already there
    }
    m_factories[type].emplace({name, []() { return static_cast<AlgorithmBase*>(new Algo()); }});
  }

  // Get a new owning pointer to an instance of an algorithm.
  // Throws if the resource isn't available, or if we aren't fully ready yet.
  template <class AlgoType> std::unique_ptr<AlgoType> get(std::string_view name) const {
    ensureReady();
    if (!available<AlgoType>(name)) {
      raise(fmt::format("No factory with name {} and type {} registered with the AlgorithmSvc",
                        name, typeid(AlgoType).name()));
    }
    const auto& factories = m_factories.at(typeid(AlgoType));
    // This creates and object and gives ownership to the unique_ptr
    return static_cast<AlgoType*>(factories.at(name)());
  }
  // Return a vector of the names of all available algorithms of a certain type
  // Just return an empty vector if no names are present (don't throw any exceptions
  // for missing entries here).
  // Throws if we aren't fully ready yet.
  template <class AlgoType> auto ls() const {
    ensureReady();
    std::vector<std::string_view> ret;
    if (available<AlgoType>()) {
      std::for_each(factory<AlgoType>().begin(), factory<AlgoType>().end(), std::back_inserter(ret),
                    [](auto&& key_value_pair) { return key_value_pair.first(); });
    }
    return ret;
  }

private:
  // Do we have any algorithms available of a certain type?
  template <class AlgoType> bool available() const { return !m_factories.count(typeid(AlgoType)); }
  // Do we have any algorithms of type AlgoType with a specific name available? Throws
  // if any of the lower-level assumptions fail
  template <class AlgoType> bool available(std::string_view name) const {
    if (name == "") {
      raise(fmt::format("Invalid name provided: '{}'", name));
    }
    return factory<AlgoType>().count(name) != 0;
  }
  template <class AlgoType> const auto& factory() const {
    if (!available<AlgoType>()) {
      raise(fmt::format("No factory for algorithm type {} provided", typeid(AlgoType).name()));
    }
    return m_factories.at(typeid(AlgoType));
  }
  // Throw an exception if we are not ready as factory to avoid giving incomplete factory
  // information to the caller
  void ensureReady() const {
    if (!ready()) {
      raise("Attempt to use AlgorithmSvc, but service not yet marked as ready");
    }
  }

  using factory_type = std::function<std::unique_ptr<AlgorithmBase>()>;

  std::map<std::type_index, std::map<std::string_view, factory_type>> m_factories;
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

