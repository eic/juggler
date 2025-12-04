# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

################################################################################
# Package: JugTrack
################################################################################

file(GLOB JugTrackPlugins_sources CONFIGURE_DEPENDS src/components/*.cpp)
gaudi_add_module(JugTrackPlugins
  SOURCES
  ${JugTrackPlugins_sources}
  LINK
  Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  EDM4HEP::edm4hep
  EDM4EIC::edm4eic
  DD4hep::DDRec
  k4FWCore::k4FWCore
  ${Acts_NAMESPACE_PREFIX}Core ${Acts_NAMESPACE_PREFIX}PluginDD4hep ${Acts_NAMESPACE_PREFIX}PluginJson
  ${ActsCore_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}ActsExamplesFramework${CMAKE_SHARED_LIBRARY_SUFFIX}
)

target_include_directories(JugTrackPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/JugBase>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_options(JugTrackPlugins PRIVATE -Wno-suggest-override)
