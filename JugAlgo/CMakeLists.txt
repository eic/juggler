# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Sylvester Joosten

################################################################################
# Package: JugAlgo
################################################################################

gaudi_add_header_only_library(JugAlgo
  LINK
  Gaudi::GaudiKernel Gaudi::GaudiAlgLib
  algorithms::core
  JugBase
)

file(GLOB JugAlgoPlugins_sources CONFIGURE_DEPENDS src/components/*.cpp)
gaudi_add_module(JugAlgoPlugins
  SOURCES
  ${JugAlgoPlugins_sources}
  LINK
  Gaudi::GaudiKernel Gaudi::GaudiAlgLib
  JugBase
  JugAlgo
  algorithms::core
)

target_include_directories(JugAlgoPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_options(JugAlgoPlugins PRIVATE -Wno-suggest-override)

install(TARGETS JugAlgo JugAlgoPlugins
  EXPORT JugAlgoTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

if(BUILD_TESTING)
  enable_testing()
endif()
