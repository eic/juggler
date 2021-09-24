cmake_minimum_required(VERSION 3.8)

project(Juggler VERSION 3.4.0)

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

include_directories(${podio_INCLUDE_DIR})
find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)
find_package(Geant4)

find_library(genfit2 genfit2 /usr/local/lib REQUIRED)

find_package(Gaudi)
add_subdirectory(JugBase)
add_subdirectory(JugDigi)
add_subdirectory(JugReco)
add_subdirectory(JugTrack)
gaudi_install(CMAKE)

# create and install Juggler.xenv file as it still has a use-case
# TODO: update workflow to not need xenv files anymore
include(cmake/xenv.cmake)
