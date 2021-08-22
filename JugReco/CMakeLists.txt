################################################################################
# Package: JugReco
################################################################################

find_package(EICD REQUIRED)
find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)

file(GLOB JugRecoPlugins_sources src/components/*.cpp)
gaudi_add_module(JugRecoPlugins
  SOURCES
  src/components/TrackingHitsCollector.cpp 
  src/components/TrackingHitsCollector2.cpp 
  src/components/TrackerHitReconstruction.cpp 
  src/components/SimpleClustering.cpp 
  src/components/CalorimeterHitReco.cpp 
  src/components/CalorimeterHitsEtaPhiProjector.cpp 
  src/components/CalorimeterHitsMerger.cpp 
  src/components/CalorimeterIslandCluster.cpp 
  src/components/CrystalEndcapsReco.cpp 
  src/components/EcalTungstenSamplingReco.cpp
  src/components/EcalTungstenSamplingCluster.cpp
  src/components/EnergyPositionClusterMerger.cpp
  src/components/ClusterRecoCoG.cpp 
  src/components/ParticlesWithTruthPID.cpp
  src/components/PhotoMultiplierReco.cpp 
  src/components/PhotoRingClusters.cpp 
  src/components/FuzzyKClusters.cpp 
  src/components/EMCalReconstruction.cpp 
  src/components/HCalReconstruction.cpp
  src/components/TopologicalCellCluster.cpp 
  src/components/ImagingPixelReco.cpp 
  src/components/ImagingPixelMerger.cpp 
  src/components/ImagingTopoCluster.cpp 
  src/components/ImagingClusterReco.cpp 
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
