#pragma once

#include <string>
#include <tuple>
#include <vector>

#include <algorithms/algorithm.h>
#include <algorithms/type_traits.h>

#include <GaudiAlg/GaudiAlgorithm.h>
#include <JugBase/DataHandle.h>

namespace Jug::Algo::detail {

// Generate properties for each of the data arguments
template <class T, bool kIsInput> class DataElement {
public:
  using value_type = std::conditional_t<algorithms::is_input_v<T>, algorithms::input_type_t<T>,
                                        algorithms::output_type_t<T>>;
  using data_type  = algorithms::data_type_t<T>;

  DataElement(gsl::not_null<GaudiAlgorithm*> owner, std::string_view name)
      : m_owner{owner}, m_data_name(m_owner, name, "") {}
  void init() {
    if (m_handle) {
      // treat error: already initialized
    }
    if (!m_data_name.empty()) {
      m_handle = std::make_unique<DataHandle<data_type>>(
          m_data_name, (kIsInput ? Gaudi::DataHandle::Reader : Gaudi::DataHandle::Writer), m_owner);
    } else if (!algorithms::is_optional_v<T>) {
      // treat error: member not optional but no collection name given
    }
  }
  value_type get() const {
    if (!m_handle) {
      return nullptr;
    }
    if constexpr (kIsInput) {
      return m_handle->get();
    } else {
      return m_handle->createAndPut();
    }
  }

private:
  GaudiAlgorithm* m_owner;
  Gaudi::Property<std::string> m_data_name;
  std::unique_ptr<DataHandle<T>> m_handle;
};

// Specialization for vectors
template <class T, class A, bool kIsInput> class DataElement<std::vector<T, A>, kIsInput> {
public:
  using value_type = std::conditional_t<algorithms::is_input_v<T>, algorithms::input_type_t<T>,
                                        algorithms::output_type_t<T>>;
  using data_type  = algorithms::data_type_t<T>;

  DataElement(gsl::not_null<GaudiAlgorithm*> owner, std::string_view name)
      : m_owner{owner}, m_data_names(m_owner, name, "") {}
  void init() {
    if (!m_handles.empty()) {
      // treat error: already initialized
    }
    if (!m_data_names.empty()) {
      for (const auto& name : m_data_names) {
        if (!name.empty()) {
          m_handles.emplace_back(std::make_unique<DataHandle<data_type>>(
              name, (kIsInput ? Gaudi::DataHandle::Reader : Gaudi::DataHandle::Writer), m_owner));
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
      if constexpr (kIsInput) {
        ret.emplace_back(handle->get());
      } else {
        ret.emplace_back(handle->createAndPut());
      }
    }
    return ret;
  }

private:
  GaudiAlgorithm* m_owner;
  Gaudi::Property<std::vector<std::string>> m_data_names;
  std::vector<std::unique_ptr<DataHandle<T>>> m_handles;
};

template <bool kIsInput, class NamesArray, class Tuple, size_t... I>
auto createElements(GaudiAlgorithm* owner, const NamesArray& names, const Tuple&,
                    std::index_sequence<I...>)
    -> std::tuple<DataElement<std::tuple_element_t<I, Tuple>, kIsInput>...> {
  return {{owner, std::get<I>(names)}...};
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
  static constexpr bool kIsInput = algorithms::is_input_v<Data>;
  using value_type               = Data;
  using data_type                = typename Data::data_type;
  constexpr static size_t kSize  = Data::kSize;
  using names_type               = typename Data::DataNames;
  using elements_type =
      decltype(createElements<kIsInput>(std::declval<GaudiAlgorithm*>(), names_type(), data_type(),
                                        std::make_index_sequence<kSize>()));

  DataProxy(gsl::not_null<GaudiAlgorithm*> owner, const names_type& names)
      : m_owner{owner}
      , m_elements{createElements<kIsInput>(m_owner, names, data_type(),
                                            std::make_index_sequence<kSize>())} {}
  void init() {
    std::apply([](auto el) { el.init(); }, m_elements);
  }
  value_type get() const {
    return getElements<value_type>(m_elements, std::make_index_sequence<kSize>());
  }

private:
  GaudiAlgorithm* m_owner;
  elements_type m_elements;
};

} // namespace Jug::Algo::detail
