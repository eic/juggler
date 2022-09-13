# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Sylvester Joosten

################################################################################
# Package: algorithms core utilities
################################################################################

set(SUBDIR "core")
set(LIBRARY "algo${SUBDIR}")
set(TARGETS ${TARGETS} ${LIBRARY} PARENT_SCOPE)

file(GLOB SRC CONFIGURE_DEPENDS src/*.cpp)

include_directories("${PROJECT_SOURCE_DIR}/${SUBDIR}/include")

add_library(${LIBRARY} SHARED ${SRC})
target_link_libraries(${LIBRARY} 
  EDM4HEP::edm4hep
  EDM4EIC::edm4eic
  DD4hep::DDRec
  fmt::fmt)
set_target_properties(${LIBRARY} PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

install(TARGETS ${LIBRARY}
  EXPORT algorithmsTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/${SUBDIR}/include/algorithms
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT dev)

# TODO: Testing
#if(BUILD_TESTING)
#  enable_testing()
#endif()

