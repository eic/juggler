cmake_minimum_required(VERSION 3.19)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

project(Juggler VERSION 4.3.0)

set(CMAKE_CXX_STANDARD 17)

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

find_library(tflite_library tensorflow-lite /opt/local/lib REQUIRED)
find_file(tflite_interpreter tensorflow/lite/interpreter.h /opt/local/include REQUIRED)
get_filename_component(tflite ${tflite_library} NAME)
get_filename_component(tflite_LIBRARY_DIR ${tflite_library} DIRECTORY)
link_directories(${tflite_LIBRARY_DIR})
get_filename_component(tflite_INCLUDE_DIR ${tflite_interpreter} DIRECTORY)
include_directories(${tflite_INCLUDE_DIR}/../..)
message(STATUS "tensorflow-lite: ${tflite} ${tflite_LIBRARY_DIR} ${tflite_INCLUDE_DIR}")

find_package(EICD REQUIRED)
find_package(NPDet REQUIRED)

find_package(podio 0.11...0.14 REQUIRED)
add_definitions("-Dpodio_VERSION_MAJOR=${podio_VERSION_MAJOR}")
add_definitions("-Dpodio_VERSION_MINOR=${podio_VERSION_MINOR}")
add_definitions("-Dpodio_VERSION_PATCH=${podio_VERSION_PATCH}")

find_package(ROOT COMPONENTS Core RIO Tree MathCore GenVector Geom REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)

find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep PluginJson)
set(Acts_VERSION_MIN "15.1.0")
set(Acts_VERSION "${Acts_VERSION_MAJOR}.${Acts_VERSION_MINOR}.${Acts_VERSION_PATCH}")
if(${Acts_VERSION} VERSION_LESS ${Acts_VERSION_MIN}
   AND NOT "${Acts_VERSION}" STREQUAL "9.9.9")
  message(FATAL_ERROR "Acts version ${Acts_VERSION_MIN} or higher required, but ${Acts_VERSION} found")
endif()
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
