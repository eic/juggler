#pragma once

#include <string>

namespace algorithms {

class PropertyBase {
public:
  PropertyBase() = default;
};

template <class TYPE> class Property : public PropertyBase {
public:
  using StorageType = TYPE;
  using ValueType   = typename std::remove_reference<StorageType>::type;

private:
  std::string m_name;
  StorageType m_value;

public:
  template <class OWNER>
  Property(OWNER* owner, const std::string& name, TYPE&& value) : m_name(name), m_value(std::forward<TYPE>(value)) {
    if (owner != nullptr) {
      owner->registerProperty(this, m_name);
    }
  }

  const ValueType& value() const { return this->m_value; }
  ValueType& value() { return const_cast<ValueType&>((const ValueType&)*this); }
};

} // namespace algorithms
