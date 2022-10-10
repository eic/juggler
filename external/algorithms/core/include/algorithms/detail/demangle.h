// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Get a demangled type name string
#pragma once

#include <cstdlib>
#include <cxxabi.h>
#include <memory>
#include <string>

#include <algorithms/error.h>

namespace algorithms::detail {

template <class T> std::string demangledName() {
  const char* mangled = typeid(T).name();
  int status          = 1; // ABI spec sets status to 0 on success, so we need a different
                           // starting value
  std::unique_ptr<char, void (*)(void*)> res{abi::__cxa_demangle(name, NULL, NULL, &status),
                                             std::free};
  if (status != 0) {
    throw algorithms::Error("Failed to demangle type name: " + mangled,
                            "algorithms::detail::demangledName()");
  }

  return res.get();
}

} // namespace algorithms::detail
