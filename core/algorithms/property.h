#pragma once

#include <any>
#include <map>
#include <string>

namespace algorithms {

// A property type that auto-registers itself with the property handler
template <class T>
class Property {
  public:
    using ValueType = T;

    // Get const reference to the value
    const ValueType& value() const { return m_value; }

    std::string_view name() const { return m_name; }

    // Only settable by explicitly calling the value operator
    // as we want the Property to mostly act as if it is constant
    template <class U> void value(U&& v) { m_value = std::forward<U>(v); }

    // act as if this is a const T
    template <typename U> bool operator==(const U& rhs) const { return m_value == m_value; }
    template <typename U> bool operator!=(const U& rhs) const { return m_value != m_value; }
    template <typename U> bool operator>(const U& rhs) const { return m_value > rhs; }
    template <typename U> bool operator>=(const U& rhs) const { return m_value >= rhs; }
    template <typename U> bool operator<(const U& rhs) const { return m_value < rhs; }
    template <typename U> bool operator<=(const U& rhs) const { return m_value <= rhs; }
    template <typename U> auto operator+(const U& rhs) const { return m_value + rhs; }
    template <typename U> auto operator-(const U& rhs) const { return m_value - rhs; }
    template <typename U> auto operator*(const U& rhs) const { return m_value * rhs; }
    template <typename U> auto operator/(const U& rhs) const { return m_value / rhs; }

  private:
    std::string_view m_name;
    T m_value;
};

// Configuration/property handling
class Configurable {
  public:

  private:
    std::map<std::string_view, std::any> m_props;
};

// Property mixin, provides all the configuration functionality for 
// our algorithms and services
// Currently an alias to Configurable
using PropertyMixin = Configurable;

}

// operator== overload not needed in C++20 as it will call the member version
#if __cpp_impl_three_way_comparison < 201711
template <class T, class U> bool operator==(const U& rhs, algorithms::Property<T>& p) {
  return p == rhs;
}
#endif
// others needed??? TODO
