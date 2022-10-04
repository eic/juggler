// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// A simple exception base class that transports both an error message and error type

#pragma once

#include <exception>
#include <string>

namespace algorithms {

class Error : public std::exception {
public:
  Error(std::string_view msg, std::string_view type = "algorithms::Error")
      : m_msg{msg}, m_type{type} {}

  virtual const char* what() const noexcept { return m_msg.c_str(); }
  virtual const char* type() const noexcept { return m_type.c_str(); }
  virtual ~Error() noexcept {}

private:
  const std::string m_msg;
  const std::string m_type;
};

} // namespace algorithms

