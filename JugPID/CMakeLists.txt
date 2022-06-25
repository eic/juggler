# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

################################################################################
# Package: JugPID
################################################################################

file(GLOB JugPIDPlugins_sources CONFIGURE_DEPENDS src/components/*.cpp)
gaudi_add_module(JugPIDPlugins
  SOURCES
  ${JugPIDPlugins_sources}
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  EDM4HEP::edm4hep
  EICD::eicd
)

target_include_directories(JugPIDPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugPIDPlugins PRIVATE -Wno-suggest-override)
