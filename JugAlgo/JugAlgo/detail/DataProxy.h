#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <algorithms/algorithm.h>
#include <algorithms/type_traits.h>

#include <GaudiAlg/GaudiAlgorithm.h>
#include <k4FWCore/DataHandle.h>

namespace Jug::Algo::detail {

enum class DataMode : unsigned { kInput, kOutput };

// Generate properties for each of the data arguments
template <class T, DataMode kMode> class DataElement {

public:
  using value_type = std::conditional_t<kMode == DataMode::kInput, algorithms::input_type_t<T>,
                                        algorithms::output_type_t<T>>;
  using data_type  = algorithms::data_type_t<T>;
  constexpr static const bool kIsOptional = algorithms::is_optional_v<T>;

  template <class Owner>
  DataElement(Owner* owner, std::string_view name)
      : m_data_name{std::make_unique<Gaudi::Property<std::string>>(owner, std::string(name), "")}
      , m_owner{owner} {}
  void init() {
    if (m_handle) {
      // treat error: already initialized
    }
    if (!m_data_name->empty()) {
      m_handle = std::make_unique<DataHandle<data_type>>(
          *m_data_name,
          ((kMode == DataMode::kInput) ? Gaudi::DataHandle::Reader : Gaudi::DataHandle::Writer),
          m_owner);
    } else if (!algorithms::is_optional_v<T>) {
      // treat error: member not optional but no collection name given
    }
  }
  value_type get() const {
    if constexpr (kIsOptional) {
      if (!m_handle) {
        return nullptr;
      }
    }
    if constexpr (kMode == DataMode::kInput) {
      return m_handle->get();
    } else {
      return m_handle->createAndPut();
    }
  }

private:
  std::unique_ptr<Gaudi::Property<std::string>>
      m_data_name; // This needs to be a pointer, else things go wrong once we go through
                   // createElements - probably something about passing the Property through an
                   // rvalue (or copy) constructor
  std::unique_ptr<DataHandle<data_type>> m_handle;
  gsl::not_null<Gaudi::Algorithm*> m_owner;
};

// Specialization for vectors
template <class T, class A, DataMode kMode> class DataElement<std::vector<T, A>, kMode> {
public:
  using value_type = std::conditional_t<kMode == DataMode::kInput, algorithms::input_type_t<T>,
                                        algorithms::output_type_t<T>>;
  using data_type  = algorithms::data_type_t<T>;

  template <class Owner>
  DataElement(Owner* owner, std::string_view name)
      : m_data_names{std::make_unique<Gaudi::Property<std::vector<std::string>>>(
            owner, std::string(name), {})}
      , m_owner{owner} {}
  void init() {
    if (!m_handles.empty()) {
      // treat error: already initialized
    }
    if (!m_data_names->empty()) {
      for (const auto& name : *m_data_names) {
        if (!name.empty()) {
          m_handles.emplace_back(std::make_unique<DataHandle<data_type>>(
              name,
              (kMode == DataMode::kInput ? Gaudi::DataHandle::Reader : Gaudi::DataHandle::Writer),
              m_owner));
        } else {
          // treat error: empty name
        }
      }
    } else {
      // OK if nothing given here, no error
    }
  }
  std::vector<value_type> get() const {
    std::vector<value_type> ret;
    for (auto& handle : m_handles) {
      if constexpr (kMode == DataMode::kInput) {
        ret.emplace_back(handle->get());
      } else {
        ret.emplace_back(handle->createAndPut());
      }
    }
    return ret;
  }

private:
  std::unique_ptr<Gaudi::Property<std::vector<std::string>>> m_data_names;
  std::vector<std::unique_ptr<DataHandle<T>>> m_handles;
  gsl::not_null<Gaudi::Algorithm*> m_owner;
};

template <DataMode kMode, class Owner, class NamesArray, class Tuple, size_t... I>
auto createElements(Owner* owner, const NamesArray& names, const Tuple&, std::index_sequence<I...>)
    -> std::tuple<DataElement<std::tuple_element_t<I, Tuple>, kMode>...> {
  return {DataElement<std::tuple_element_t<I, Tuple>, kMode>(owner, std::get<I>(names))...};
}

// Call ::get() on each element of the HandleTuple, and return the result in the format of
// ReturnTuple
template <class ReturnTuple, class HandleTuple, size_t... I>
ReturnTuple getElements(HandleTuple& handles, std::index_sequence<I...>) {
  return {std::get<I>(handles).get()...};
}

// Create data handle structure for all members

template <class Data> class DataProxy {
public:
  static constexpr DataMode kMode =
      (algorithms::is_input_v<Data> ? DataMode::kInput : DataMode::kOutput);
  using value_type              = typename Data::value_type;
  using data_type               = typename Data::data_type;
  constexpr static size_t kSize = Data::kSize;
  using names_type              = typename Data::key_type;
  using elements_type =
      decltype(createElements<kMode>(std::declval<GaudiAlgorithm*>(), names_type(), data_type(),
                                     std::make_index_sequence<kSize>()));

  template <class Owner>
  DataProxy(Owner* owner, const names_type& names)
      : m_elements{
            createElements<kMode>(owner, names, data_type(), std::make_index_sequence<kSize>())} {}
  void init() {
    std::apply([](auto&&... el) { (el.init(), ...); }, m_elements);
  }
  value_type get() const {
    return getElements<value_type>(m_elements, std::make_index_sequence<kSize>());
  }

private:
  elements_type m_elements;
};

} // namespace Jug::Algo::detail
