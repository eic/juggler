cmake_minimum_required(VERSION 3.19)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

project(Juggler VERSION 4.1.0)

set(CMAKE_CXX_STANDARD 17)

find_package(EICD REQUIRED)
find_package(NPDet REQUIRED)

find_package(podio 0.11...0.14 REQUIRED)
add_definitions("-Dpodio_VERSION_MAJOR=${podio_VERSION_MAJOR}")
add_definitions("-Dpodio_VERSION_MINOR=${podio_VERSION_MINOR}")
add_definitions("-Dpodio_VERSION_PATCH=${podio_VERSION_PATCH}")

find_package(ROOT COMPONENTS Core RIO Tree MathCore GenVector Geom REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)

find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep PluginJson)
add_definitions("-DActs_VERSION_MAJOR=${Acts_VERSION_MAJOR}")
add_definitions("-DActs_VERSION_MINOR=${Acts_VERSION_MINOR}")
add_definitions("-DActs_VERSION_PATCH=${Acts_VERSION_PATCH}")

find_library(genfit2 genfit2 /usr/local/lib REQUIRED)
find_path(genfit2_INCLUDE_DIR NAMES GFGbl.h PATHS /usr/local/include ${genfit2}/../include REQUIRED)

find_package(Gaudi)
add_subdirectory(JugBase)
add_subdirectory(JugDigi)
add_subdirectory(JugFast)
add_subdirectory(JugPID)
add_subdirectory(JugReco)
add_subdirectory(JugTrack)
gaudi_install(CMAKE)

# create and install Juggler.xenv file as it still has a use-case
# TODO: update workflow to not need xenv files anymore
include(cmake/xenv.cmake)
