################################################################################
# Package: JugDigi
################################################################################

file(GLOB JugDigiPlugins_sources src/components/*.cpp)
gaudi_add_module(JugDigiPlugins
  SOURCES
  src/components/CalorimeterHitDigi.cpp 
  src/components/CalorimeterBirksCorr.cpp
  src/components/PhotoMultiplierDigi.cpp 
  src/components/UFSDTrackerDigi.cpp 
  src/components/SiliconTrackerDigi.cpp 
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  EDM4HEP::edm4hep
  EICD::eicd
)

target_include_directories(JugDigiPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugDigiPlugins PRIVATE -Wno-suggest-override)
