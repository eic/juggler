# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Wouter Deconinck, Whitney Armstrong

################################################################################
# Package: JugBase
################################################################################

gaudi_add_header_only_library(JugBase
  LINK
  Gaudi::GaudiKernel
  podio::podioRootIO
  ROOT::Core ROOT::RIO ROOT::Tree
  DD4hep::DDRec
  ActsCore ActsPluginDD4hep
  k4FWCore::k4FWCore
)

file(GLOB JugBasePlugins_sources CONFIGURE_DEPENDS src/components/*.cpp)
gaudi_add_module(JugBasePlugins
  SOURCES
  ${JugBasePlugins_sources}
  LINK
  Gaudi::GaudiKernel
  ROOT::Core ROOT::RIO ROOT::Tree
  JugBase
  EDM4HEP::edm4hep
  DD4hep::DDRec
  ActsCore ActsPluginDD4hep ActsPluginJson
  EDM4EIC::edm4eic
)

target_include_directories(JugBasePlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

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
