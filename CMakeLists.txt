# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Wouter Deconinck, Whitney Armstrong

cmake_minimum_required(VERSION 3.19)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

project(Juggler VERSION 4.3.0)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
if(NOT CMAKE_CXX_STANDARD MATCHES "17|20")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands as json for run-clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Also use clang-tidy integration in CMake
option(ENABLE_CLANG_TIDY "Enable clang-tidy integration in cmake" OFF)
if(ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
  if (CLANG_TIDY_EXE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}" CACHE STRING "" FORCE)
  else()
    set(CMAKE_CXX_CLANG_TIDY "" CACHE STRING "" FORCE)
  endif()
endif()

# Set default build type
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "RelWithDebInfo")
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Set all warnings
if(NOT CMAKE_BUILD_TYPE MATCHES Release)
  add_compile_options(-Wall -Wextra -Werror)
endif()

find_package(Microsoft.GSL CONFIG)

find_package(EDM4EIC REQUIRED)
find_package(EDM4HEP 0.4.1 REQUIRED)

find_package(podio 0.14.1 REQUIRED)
add_definitions("-Dpodio_VERSION_MAJOR=${podio_VERSION_MAJOR}")
add_definitions("-Dpodio_VERSION_MINOR=${podio_VERSION_MINOR}")
add_definitions("-Dpodio_VERSION_PATCH=${podio_VERSION_PATCH}")

find_package(ROOT COMPONENTS Core RIO Tree MathCore GenVector Geom REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)

find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep PluginJson)
set(Acts_VERSION_MIN "20.2.0")
set(Acts_VERSION "${Acts_VERSION_MAJOR}.${Acts_VERSION_MINOR}.${Acts_VERSION_PATCH}")
if(${Acts_VERSION} VERSION_LESS ${Acts_VERSION_MIN}
   AND NOT "${Acts_VERSION}" STREQUAL "9.9.9")
  message(FATAL_ERROR "Acts version ${Acts_VERSION_MIN} or higher required, but ${Acts_VERSION} found")
endif()
add_definitions("-DActs_VERSION_MAJOR=${Acts_VERSION_MAJOR}")
add_definitions("-DActs_VERSION_MINOR=${Acts_VERSION_MINOR}")
add_definitions("-DActs_VERSION_PATCH=${Acts_VERSION_PATCH}")

find_package(IRT REQUIRED)

## algorithms dependency
option(USE_SYSTEM_ALGORITHMS "Use system-provided algorithms" OFF)
if (USE_SYSTEM_ALGORITHMS)
  find_package(algorithms)
else()
  add_subdirectory(external/algorithms)
  add_library(algorithms::core ALIAS algocore)
  add_library(algorithms::truth ALIAS algotruth)
  add_library(algorithms::calorimetry ALIAS algocalorimetry)
endif()

find_package(Gaudi)
add_subdirectory(JugBase)
add_subdirectory(JugAlgo)
add_subdirectory(JugDigi)
add_subdirectory(JugFast)
add_subdirectory(JugPID)
add_subdirectory(JugReco)
add_subdirectory(JugTrack)
gaudi_install(CMAKE)


# create and install Juggler.xenv file as it still has a use-case
# TODO: update workflow to not need xenv files anymore
include(cmake/xenv.cmake)
