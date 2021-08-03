################################################################################
# Package: JugTrack
################################################################################

find_package(EICD REQUIRED)
find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep )
find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)
find_library(genfit2 genfit2 /usr/local/lib REQUIRED)

file(GLOB JugTrackPlugins_sources src/components/*.cpp)
gaudi_add_module(JugTrackPlugins
  SOURCES
  src/components/GenFitTrackFitter.cpp 
  src/components/SingleTrackSourceLinker.cpp 
  src/components/TrackerSourceLinker.cpp 
  src/components/Tracker2SourceLinker.cpp 
  src/components/TrackerSourcesLinker.cpp 
  src/components/TrackFindingAlgorithm.cpp 
  src/components/TrackFindingAlgorithmFunction.cpp 
  src/components/TrackFittingAlgorithm.cpp 
  src/components/TrackFittingFunction.cpp 
  src/components/TestACTSLogger.cpp 
  src/components/TrackParamTruthInit.cpp 
  src/components/TrackParamClusterInit.cpp 
  src/components/TrackParamVertexClusterInit.cpp 
  src/components/TrackParamImagingClusterInit.cpp 
  src/components/ParticlesFromTrackFit.cpp 
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  NPDet::DD4podIO
  EICD::eicd
  DD4hep::DDRec
  ActsCore
  ${genfit2}
)

target_include_directories(JugTrackPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/JugBase>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

target_compile_options(JugTrackPlugins PRIVATE -Wno-suggest-override)

#gaudi_add_test(loggerTest
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/testLogger.py
#               )
