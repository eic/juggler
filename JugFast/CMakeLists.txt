# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

################################################################################
# Package: JugFast
################################################################################

gaudi_add_module(JugFastPlugins
  SOURCES
  src/components/InclusiveKinematicsTruth.cpp
  src/components/ClusterMerger.cpp 
  src/components/MatchClusters.cpp 
  src/components/MC2SmearedParticle.cpp 
  src/components/ParticlesWithTruthPID.cpp 
  src/components/SmearedFarForwardParticles.cpp 
  src/components/TruthEnergyPositionClusterMerger.cpp
  src/components/TruthClustering.cpp
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  EDM4HEP::edm4hep
  EICD::eicd
)

target_include_directories(JugFastPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugFastPlugins PRIVATE -Wno-suggest-override)
