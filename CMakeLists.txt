cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 17)

find_package(EICD REQUIRED)
#set(CMAKE_PREFIX_PATH $ENV{HOME}/lib CMAKE_PREFIX_PATH)
find_package(NPDet  REQUIRED)
#---------------------------------------------------------------
# to use ROOT targets, find_package ROOT must come before find_package Gaudi
# (for now) # see issue Gaudi#103
set(PODIO $ENV{PODIO})
set(CMAKE_MODULE_PATH CMAKE_MODULE_PATH PODIO)
find_package(podio 0.11.0 REQUIRED)

## There was a breaking change in how collections are stored in podio between
## version v0.13.0 and v0.13.1. This checks for podio version and sets the 
## PODIO_BEFORE_0_13_1 preprocessor variable
if (${podio_VERSION_MAJOR} EQUAL 0 AND ${podio_VERSION_MINOR} LESS_EQUAL 13)
  if (${podio_VERSION_MINOR} EQUAL 13 AND ${podio_VERSION_PATCH} EQUAL 0)
    add_definitions("-DPODIO_BEFORE_0_13_1")
  elseif (${podio_VERSION_MINOR} LESS 13)
    add_definitions("-DPODIO_BEFORE_0_13_1")
  endif()
endif()

include_directories(${podio_INCLUDE_DIR})
find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)
find_package(Geant4)

find_library(genfit2 genfit2 /usr/local/lib REQUIRED)

find_package(GaudiProject)

gaudi_project(Juggler v1r0
              USE Gaudi v34r0)

install(FILES ${CMAKE_BINARY_DIR}/run 
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ  WORLD_EXECUTE 
  DESTINATION bin/
  )
