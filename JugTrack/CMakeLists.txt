################################################################################
# Package: JugTrack
################################################################################
gaudi_subdir(JugTrack v1r0)

find_package(EICD REQUIRED)
#set(PODIO $ENV{PODIO})
#set(CMAKE_MODULE_PATH CMAKE_MODULE_PATH PODIO)
#find_package(podio 0.11.01 REQUIRED)
#include_directories(${podio_INCLUDE_DIR})

find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep )
find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)

# this declaration will not be needed in the future
gaudi_depends_on_subdirs(GaudiAlg GaudiKernel)

#gaudi_add_library(JugBase
#  src/*.cpp 
#                  INCLUDE_DIRS EICD PODIO ROOT $ENV{HOME}/stow/podio/include
#                  LINK_LIBRARIES GaudiAlgLib GaudiKernel ROOT DD4hep::DDG4IO
#                  PUBLIC_HEADERS JugBase)
#target_link_libraries(JugBase 
#  podio::podioRootIO
#  )

gaudi_add_module(JugTrackPlugins
  src/components/TrackerSourceLinker.cpp 
  src/components/Tracker2SourceLinker.cpp 
  src/components/TrackerSourcesLinker.cpp 
  src/components/TrackFindingAlgorithm.cpp 
  src/components/TrackFindingAlgorithmFunction.cpp 
  #src/components/TrackFittingAlgorithm.cpp 
  #src/components/TrackFittingFunction.cpp 
  src/components/TestACTSLogger.cpp 
  src/components/TrackParamTruthInit.cpp 
  src/components/TrackParamClusterInit.cpp 
  src/components/TrackParamVertexClusterInit.cpp 
  src/components/ParticlesFromTrackFit.cpp 
  LINK_LIBRARIES GaudiAlgLib GaudiKernel JugBase ROOT NPDet::DD4podIO EICD::eicd DDRec Acts)

target_compile_options(JugTrackPlugins PRIVATE -Wno-suggest-override)

#gaudi_add_test(loggerTest
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/testLogger.py
#               )
#gaudi_add_test(CheckReadCollectionSize
#               ENVIRONMENT PYTHONPATH+=${PODIO_PYTHON_DIR}
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               COMMAND python FWCore/tests/scripts/check_coll_after_read.py
#               DEPENDS ReadTest)