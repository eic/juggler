################################################################################
# Package: JugPID
################################################################################

include_directories( 
  # For now assume that newly installed eicd::(...) event structures are available 
  # in the installation area;
  ${CMAKE_INSTALL_PREFIX}/include 

  ${PROJECT_SOURCE_DIR}/JugTrack
)

gaudi_add_module(JugPIDPlugins
  SOURCES
  src/components/TOFPID.cpp 
  src/components/FuzzyKClusters.cpp 
  src/components/PhotoRingClusters.cpp
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  NPDet::DD4podIO
  EICD::eicd
)

target_include_directories(JugPIDPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugPIDPlugins PRIVATE -Wno-suggest-override)
