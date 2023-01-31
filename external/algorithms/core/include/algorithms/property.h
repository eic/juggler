// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Defines the PropertyMixin
// This base class provide access to the Property<T> object, a self-registering
// configurable property that acts as a bare object T from a performance point-of-view
//
#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <variant>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithms/detail/upcast.h>
#include <algorithms/error.h>
#include <algorithms/name.h>

namespace algorithms {

class PropertyError : public Error {
public:
  PropertyError(std::string_view msg) : Error{msg, "algorithms::PropertyError"} {}
};

// Data types supported for Properties, defined as std::variant. This allows for
// automatic Property registration with the calling framework, and enables compile-time
// type errors.
using PropertyValue =
    std::variant<bool, uint32_t, int32_t, uint64_t, int64_t, double, std::string, std::vector<bool>,
                 std::vector<uint32_t>, std::vector<int32_t>, std::vector<uint64_t>,
                 std::vector<int64_t>, std::vector<double>, std::vector<std::string>>;

// Configuration/property handling
class PropertyMixin {
public:
  class PropertyHandle;
  using PropertyMap = std::map<std::string_view, PropertyHandle&>;

  template <typename T> void setProperty(std::string_view name, T&& value) {
    m_props.at(name).set(static_cast<detail::upcast_type_t<T>>(value));
  }
  template <typename T> T getProperty(std::string_view name) const {
    return std::get<T>(m_props.at(name).get());
  }
  const PropertyMap& getProperties() const { return m_props; }
  bool hasProperty(std::string_view name) const {
    return m_props.count(name) && m_props.at(name).hasValue();
  }
  // get a vector of the names of all missing (unset) properties
  auto missingProperties() const {
    std::vector<std::string_view> missing;
    for (const auto& [name, prop] : m_props) {
      if (!prop.hasValue()) {
        missing.push_back(name);
      }
    }
    return missing;
  }
  // Throw an exception if any properties are not set
  void validate() const {
    const auto missing = missingProperties();
    if (!missing.empty()) {
      throw PropertyError(fmt::format("Missing properties: {}", missing));
    }
  }

private:
  void registerProperty(PropertyHandle& prop) {
    if (m_props.count(prop.name())) {
      throw PropertyError(fmt::format("Duplicate property name: {}", prop.name()));
    }
    m_props.emplace(prop.name(), prop);
  }

  PropertyMap m_props;

public:
  class PropertyHandle : public NameMixin {
  public:
    PropertyHandle(std::string_view name, std::string_view description)
        : NameMixin{name, description} {}
    virtual void set(const PropertyValue& v) = 0;
    virtual PropertyValue get() const        = 0;
    bool hasValue() const { return m_has_value; }

  protected:
    bool m_has_value = false;
  };

  // A property type that auto-registers itself with the property handler
  // Essentially a simplified and const-like version of Gaudi::Property
  template <class T> class Property : public PropertyHandle {
  public:
    using value_type = T;
    using impl_type  = detail::upcast_type_t<T>;

    Property(PropertyMixin* owner, std::string_view name, std::string_view description)
        : PropertyHandle{name, description} {
      if (owner) {
        owner->registerProperty(*this);
      } else {
        throw PropertyError(
            fmt::format("Attempting to create Property '{}' without valid owner", name));
      }
    }
    Property(PropertyMixin* owner, std::string_view name, const value_type& v,
             std::string_view description)
        : Property(owner, name, description) {
      set(static_cast<impl_type>(v));
    }

    Property()                = delete;
    Property(const Property&) = default;
    Property& operator=(const Property&) = default;

    // Only settable by explicitly calling the ::set() member function
    // as we want the Property to mostly act as if it is constant
    virtual void set(const PropertyValue& v) {
      m_value     = static_cast<value_type>(std::get<impl_type>(v));
      m_has_value = true;
    }
    // virtual getter for use from PropertyHandle - use ::value() instead for direct member
    // access
    virtual PropertyValue get() const { return static_cast<impl_type>(m_value); }

    // Direct access to the value. Use this one whenever possible (or go through the
    // automatic casting)
    const value_type& value() const { return m_value; }

    // automatically cast to T
    operator T() const { return m_value; }

    // act as if this is a const T
    template <typename U> bool operator==(const U& rhs) const { return m_value == rhs; }
    template <typename U> bool operator!=(const U& rhs) const { return m_value != rhs; }
    template <typename U> bool operator>(const U& rhs) const { return m_value > rhs; }
    template <typename U> bool operator>=(const U& rhs) const { return m_value >= rhs; }
    template <typename U> bool operator<(const U& rhs) const { return m_value < rhs; }
    template <typename U> bool operator<=(const U& rhs) const { return m_value <= rhs; }
    template <typename U> decltype(auto) operator+(const U& rhs) const { return m_value + rhs; }
    template <typename U> decltype(auto) operator-(const U& rhs) const { return m_value - rhs; }
    template <typename U> decltype(auto) operator*(const U& rhs) const { return m_value * rhs; }
    template <typename U> decltype(auto) operator/(const U& rhs) const { return m_value / rhs; }

    // stl collection helpers if needed
    // forced to be templated so we only declare them when used
    template <class U = const value_type> decltype(auto) size() const { return value().size(); }
    template <class U = const value_type> decltype(auto) length() const { return value().length(); }
    template <class U = const value_type> decltype(auto) empty() const { return value().empty(); }
    template <class U = value_type> decltype(auto) clear() { value().clear(); }
    template <class U = const value_type> decltype(auto) begin() const { return value().begin(); }
    template <class U = const value_type> decltype(auto) end() const { return value().end(); }
    template <class U = value_type> decltype(auto) begin() { return value().begin(); }
    template <class U = value_type> decltype(auto) end() { return value().end(); }
    template <class Arg> decltype(auto) operator[](const Arg& arg) const { return value()[arg]; }
    template <class Arg> decltype(auto) operator[](const Arg& arg) { return value()[arg]; }
    template <class U = const value_type>
    decltype(auto) find(const typename U::key_type& key) const {
      return value().find(key);
    }
    template <class U = value_type> decltype(auto) find(const typename U::key_type& key) {
      return value().find(key);
    }
    template <class Arg> decltype(auto) erase(const Arg& arg) { return value().erase(arg); }

    // In case our property has operator (), delegate the operator()
    template <class... Args>
    decltype(std::declval<value_type>()(std::declval<Args&&>()...))
    operator()(Args&&... args) const {
      return m_value()(std::forward<Args>(args)...);
    }

  private:
    T m_value;
  };
}; // namespace algorithms

} // namespace algorithms

// operator== overload not needed in C++20 as it will call the member version
#if __cpp_impl_three_way_comparison < 201711
template <class T, class U>
bool operator==(const U& rhs, const algorithms::PropertyMixin::Property<T>& p) {
  return p == rhs;
}
#endif
template <class T>
std::ostream& operator<<(std::ostream& os, const algorithms::PropertyMixin::Property<T>& p) {
  return os << p.value();
}

// Make Property formateble
template <class T>
struct fmt::formatter<algorithms::PropertyMixin::Property<T>> : fmt::formatter<T> {};
// others needed??? TODO
