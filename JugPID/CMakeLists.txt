# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

################################################################################
# Package: JugPID
################################################################################

file(GLOB JugPIDPlugins_sources src/components/*.cpp)

# TODO: is DD4hep::DDRec needed (added for debugging missing algorithm issue)
gaudi_add_module(JugPIDPlugins
  SOURCES
  src/components/FuzzyKClusters.cpp 
  src/components/PhotoRingClusters.cpp
  src/components/TestIRTAlgorithm.cpp
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  EDM4HEP::edm4hep
  EICD::eicd
  DD4hep::DDRec
  IRT
)

# TODO: is ${CMAKE_SOURCE_DIR}/JugBase needed? (added for debugging missing algo issue)
# TODO: need to include also ${IRT_INCLUDE_DIR}>/IRT, since `#include <...>` calls in IRT
#       headers do not have `IRT/` in their path
target_include_directories(JugPIDPlugins PUBLIC
  $<BUILD_INTERFACE:${IRT_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${IRT_INCLUDE_DIR}>/IRT
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/JugBase>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugPIDPlugins PRIVATE -Wno-suggest-override)
