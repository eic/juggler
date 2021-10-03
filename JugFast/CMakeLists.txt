################################################################################
# Package: JugFast
################################################################################

find_package(EICD REQUIRED)

find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)

gaudi_add_module(JugFastPlugins
  SOURCES
  src/components/MatchClusters.cpp 
  src/components/MC2SmearedParticle.cpp 
  src/components/ParticlesWithTruthPID.cpp 
  src/components/SmearedFarForwardParticles.cpp 
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  NPDet::DD4podIO
  EICD::eicd
)

target_include_directories(JugFastPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugFastPlugins PRIVATE -Wno-suggest-override)
