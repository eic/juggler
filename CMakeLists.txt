cmake_minimum_required(VERSION 3.15)
project(Juggler 
  LANGUAGES CXX
  VERSION 0.5.0
  DESCRIPTION "Juggler of reconstruction"
  HOMEPAGE_URL "https://eicweb.phy.anl.gov/EIC/juggler"
  )

set(CMAKE_PREFIX_PATH $ENV{HOME}/lib CMAKE_PREFIX_PATH)
find_package(NPDet  REQUIRED)
#---------------------------------------------------------------
# to use ROOT targets, find_package ROOT must come before find_package Gaudi
# (for now) # see issue Gaudi#103
set(PODIO $ENV{PODIO})
set(CMAKE_MODULE_PATH CMAKE_MODULE_PATH PODIO)
find_package(podio 0.11.0 REQUIRED)
include_directories(${podio_INCLUDE_DIR})
find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)
find_package(Geant4)


find_package(Gaudi)

add_subdirectory(JugBase)
add_subdirectory(JugDigi)
add_subdirectory(JugReco)



install(FILES ${CMAKE_BINARY_DIR}/run 
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ  WORLD_EXECUTE 
  DESTINATION bin/
  )
