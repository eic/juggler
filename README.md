Project Juggler
===============

Concurrent event processor for NP experiments, based on the Gaudi framework.

Overview
--------

### Running Juggler

Here is an example for topside
```
./scripts/run_topside.py -i inputevents.hepmc -o output_events.root -n 10
../where_ever/../juggler/build/run gaudirun.py options/example_reconstruction.py
```

### `JugBase`

### `JugDigi`

### `JugReco`

#### Hit reconstruction algorithms

#### Source Linkers

#### Initial Track Parameters 

#### Track Fitting

```
CalorimeterIslandCluster.cpp
ClusterRecoCoG.cpp
CrystalEndcapsReco.cpp
EMCalReconstruction.cpp
FuzzyKClusters.cpp
ParticlesFromTrackFit.cpp
PhotoMultiplierReco.cpp
PhotoRingClusters.cpp
SimpleClustering.cpp
TestACTSLogger.cpp
TrackerHitReconstruction.cpp
TrackerSourceLinker.cpp
TrackFindingAlgorithm.cpp
TrackingHitsSourceLinker.cpp
TrackParamClusterInit.cpp
TrackParamTruthInit.cpp
TrackParamVertexClusterInit.cpp
```
