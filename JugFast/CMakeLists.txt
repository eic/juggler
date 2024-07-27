# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

################################################################################
# Package: JugFast
################################################################################

file(GLOB JugFastPlugins_sources CONFIGURE_DEPENDS src/components/*.cpp)
gaudi_add_module(JugFastPlugins
  SOURCES
  ${JugFastPlugins_sources}
  LINK
  Gaudi::GaudiKernel
  JugBase
  JugAlgo
  EICrecon::algorithms_calorimetry_library
  EICrecon::algorithms_reco_library
  algorithms::core algorithms::truth
  ROOT::Core ROOT::RIO ROOT::Tree
  EDM4HEP::edm4hep
  EDM4EIC::edm4eic
  k4FWCore::k4FWCore
  GENCONF_OPTIONS --load-library=JugAlgoPreload.so
)

target_include_directories(JugFastPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugFastPlugins PRIVATE -Wno-suggest-override)
