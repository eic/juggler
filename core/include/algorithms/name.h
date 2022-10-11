// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Defines NameMixin - simple base class to provide name functionality
//
#pragma once

#include <string>

namespace algorithms {

// Simple name Mixin providing consistent name API
class NameMixin {
public:
  NameMixin(std::string_view name, std::string_view description)
      : m_name{name}, m_description{description} {}
  std::string_view name() const { return m_name; }
  std::string_view description() const { return m_description; }

private:
  const std::string m_name;
  const std::string m_description;
};

} // namespace algorithms

