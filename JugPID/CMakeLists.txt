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
  NPDet::DD4podIO
  EICD::eicd
  DD4hep::DDRec
)

# TODO: is ${CMAKE_SOURCE_DIR}/JugBase needed? (added for debugging missing algo issue)
target_include_directories(JugPIDPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/JugBase>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugPIDPlugins PRIVATE -Wno-suggest-override)
