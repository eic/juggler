################################################################################
# Package: JugDigi
################################################################################

find_package(EICD REQUIRED)
#set(PODIO $ENV{PODIO})
#set(CMAKE_MODULE_PATH CMAKE_MODULE_PATH PODIO)
#find_package(podio 0.11.01 REQUIRED)
#include_directories(${podio_INCLUDE_DIR})

find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)

file(GLOB JugDigiPlugins_sources src/components/*.cpp)
gaudi_add_module(JugDigiPlugins
  SOURCES
  src/components/CalorimeterHitDigi.cpp 
  src/components/CrystalEndcapsDigi.cpp 
  src/components/EMCalorimeterDigi.cpp 
  src/components/EcalTungstenSamplingDigi.cpp 
  src/components/ExampleCaloDigi.cpp 
  src/components/HadronicCalDigi.cpp 
  src/components/PhotoMultiplierDigi.cpp 
  src/components/SOIPIXTrackerDigi.cpp 
  src/components/UFSDTrackerDigi.cpp 
  src/components/SiliconTrackerDigi.cpp 
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  NPDet::DD4podIO
  EICD::eicd
)

target_include_directories(JugDigiPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugDigiPlugins PRIVATE -Wno-suggest-override)
