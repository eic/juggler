# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Wouter Deconinck, Whitney Armstrong, Chao Peng, Sylvester Joosten, Chao

################################################################################
# Package: JugReco
################################################################################

file(GLOB JugRecoPlugins_sources CONFIGURE_DEPENDS src/components/*.cpp)
gaudi_add_module(JugRecoPlugins
  SOURCES
  ${JugRecoPlugins_sources}
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase JugAlgo
  algorithms::core algorithms::calorimetry
  ROOT::Core ROOT::RIO ROOT::Tree
  EDM4HEP::edm4hep
  EDM4EIC::edm4eic
  DD4hep::DDRec
)

target_include_directories(JugRecoPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugRecoPlugins PRIVATE -Wno-suggest-override)
