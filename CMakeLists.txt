cmake_minimum_required(VERSION 3.2.0)

#---------------------------------------------------------------
# to use ROOT targets, find_package ROOT must come before find_package Gaudi
# (for now) # see issue Gaudi#103
set(PODIO $ENV{PODIO})
set(CMAKE_MODULE_PATH CMAKE_MODULE_PATH PODIO)
find_package(podio 0.11.01 REQUIRED)
include_directories(${podio_INCLUDE_DIR})
find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)
find_package(Geant4)

find_package(GaudiProject)

gaudi_project(Juggler v1r0
              USE Gaudi v33r1)
