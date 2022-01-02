################################################################################
# Package: JugReco
################################################################################

file(GLOB JugRecoPlugins_sources src/components/*.cpp)
gaudi_add_module(JugRecoPlugins
  SOURCES
  src/components/FarForwardParticles.cpp
  src/components/FarForwardParticlesOMD.cpp
  src/components/TrackingHitsCollector.cpp 
  src/components/TrackingHitsCollector2.cpp 
  src/components/TrackerHitReconstruction.cpp 
  src/components/SimpleClustering.cpp 
  src/components/CalorimeterHitReco.cpp 
  src/components/CalorimeterHitsEtaPhiProjector.cpp 
  src/components/CalorimeterHitsMerger.cpp 
  src/components/CalorimeterIslandCluster.cpp 
  src/components/EnergyPositionClusterMerger.cpp
  src/components/ClusterRecoCoG.cpp 
  src/components/ParticleCollector.cpp
  src/components/ImagingPixelReco.cpp 
  src/components/ImagingPixelMerger.cpp 
  src/components/ImagingPixelDataCombiner.cpp 
  src/components/ImagingPixelDataSorter.cpp 
  src/components/ImagingTopoCluster.cpp 
  src/components/ImagingClusterReco.cpp 
  src/components/InclusiveKinematicsElectron.cpp
  src/components/InclusiveKinematicsJB.cpp
  src/components/InclusiveKinematicsDA.cpp
  src/components/InclusiveKinematicsSigma.cpp
  src/components/InclusiveKinematicseSigma.cpp
  src/components/PhotoMultiplierReco.cpp
  LINK
  Gaudi::GaudiAlgLib Gaudi::GaudiKernel
  JugBase
  ROOT::Core ROOT::RIO ROOT::Tree
  NPDet::DD4podIO
  EICD::eicd
  DD4hep::DDRec
)

target_include_directories(JugRecoPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_options(JugRecoPlugins PRIVATE -Wno-suggest-override)
