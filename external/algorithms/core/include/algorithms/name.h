// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Defines NameMixin - simple base class to provide name functionality
//
#pragma once

#include <cstdint>
#include <functional> // for std::hash
#include <limits>
#include <string>

namespace algorithms {

// Simple name Mixin providing consistent name API
class NameMixin {
public:
  using IdType = uint32_t;

  NameMixin(std::string_view name, std::string_view description)
      : m_name{name}, m_description{description}, m_id{calc_hash(name)} {}
  std::string_view name() const { return m_name; }
  std::string_view description() const { return m_description; }
  IdType id() const { return m_id; }

private:
  // Helper function to initialize the unique ID based on this entity's name 
  IdType calc_hash(std::string_view s) {
    std::hash<std::string_view> hash_alg;
    const auto fullID = hash_alg(s);
    const IdType max  = std::numeric_limits<IdType>::max();
    return static_cast<IdType>(fullID & max);
  }

  const std::string m_name;
  const std::string m_description;
  const uint32_t m_id;
};

} // namespace algorithms

