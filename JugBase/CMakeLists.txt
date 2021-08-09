################################################################################
# Package: JugBase
################################################################################

find_package(EICD)
find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep )
find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)
find_library(genfit2 genfit2 /usr/local/lib REQUIRED)

gaudi_add_library(JugBase
  SOURCES
  src/ACTSLogger.cpp
  src/PodioDataSvc.cpp
  src/KeepDropSwitch.cpp
  src/Utilities/Helpers.cpp
  src/Utilities/Paths.cpp
  src/Utilities/Options.cpp
  src/Plugins/DD4hepBField.cpp
  LINK
  Gaudi::GaudiKernel Gaudi::GaudiAlgLib
  podio::podioRootIO
  ROOT::Core ROOT::RIO ROOT::Tree
  DD4hep::DDG4IO DD4hep::DDRec
  ActsCore ActsPluginDD4hep
  ${genfit2}
)

target_include_directories(JugBase PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugBase PRIVATE -Wno-suggest-override)

file(GLOB JugBasePlugins_sources src/components/*.cpp)
gaudi_add_module(JugBasePlugins
  SOURCES
  src/components/EICDataSvc.cpp
  src/components/GeoSvc.cpp
  src/components/ParticleSvc.cpp
  src/components/InputCopier.cpp
  src/components/MC2DummyParticle.cpp
  src/components/PodioInput.cpp
  src/components/PodioOutput.cpp
  src/components/type.cpp
  src/components/ReadTestConsumer.cxx
  LINK
  Gaudi::GaudiKernel Gaudi::GaudiAlgLib
  GaudiKernel GaudiAlgLib
  ROOT::Core ROOT::RIO ROOT::Tree
  JugBase
  NPDet::DD4podIO
  DD4hep::DDRec
  ActsCore ActsPluginDD4hep
  EICD::eicd
  ${genfit2}
)

target_include_directories(JugBasePlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugBasePlugins PRIVATE -Wno-suggest-override)

install(TARGETS JugBase JugBasePlugins
  EXPORT JugBaseTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

if(BUILD_TESTING)
  enable_testing()
endif()

#add_test(NAME ProduceForReadTest
#         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#         COMMAND ${CMAKE_BINARY_DIR}/run ${PROJECT_SOURCE_DIR}/JugBase/scripts/gaudirun tests/options/simple_producer.py)
#add_test(NAME ReadTest
#         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#         COMMAND ${CMAKE_BINARY_DIR}/run ${PROJECT_SOURCE_DIR}/JugBase/scripts/gaudirun tests/options/simple_reader.py)
#add_test(NAME ReadGeoTest
#         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#         COMMAND ${CMAKE_BINARY_DIR}/run ${PROJECT_SOURCE_DIR}/JugBase/scripts/gaudirun tests/options/reader_with_geosvc.py)
#add_test(NAME CheckReadCollectionSize
#         ENVIRONMENT PYTHONPATH+=${PODIO_PYTHON_DIR}
#         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#         COMMAND python JugBase/tests/scripts/check_coll_after_read.py
#         DEPENDS ReadTest)
