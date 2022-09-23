// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Defines NameMixin - simple base class to provide name functionality
//
#pragma once

namespace algorithms {

// Simple name Mixin providing consistent name API
class NameMixin {
public:
  NameMixin(std::string_view name) : m_name{name} {}
  std::string_view name() const { return m_name; }

private:
  const std::string m_name;
};

} // namespace algorithms

