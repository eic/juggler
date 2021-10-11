################################################################################
# Package: JugTrack
################################################################################

file(GLOB JugTrackPlugins_sources src/components/*.cpp)
gaudi_add_module(JugTrackPlugins
  SOURCES
  src/components/ConformalXYPeakProtoTracks.cpp
  src/components/FinderAlgoTemplate.cpp 
  src/components/ProtoTrackMatching.cpp 
  src/components/GenFitTrackFitter.cpp 
  src/components/SingleTrackSourceLinker.cpp 
  src/components/TrackerSourceLinker.cpp 
  src/components/CKFTracking.cpp 
  src/components/CKFTrackingFunction.cpp 
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
  src/components/ProjectTrack.cpp 
  src/components/TrajectoryFromTrackFit.cpp 
  src/components/TruthTrackSeeding.cpp 
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
  ${genfit2_INCLUDE_DIR}
)

target_compile_options(JugTrackPlugins PRIVATE -Wno-suggest-override)


#gaudi_add_test(loggerTest
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/testLogger.py
#               )
