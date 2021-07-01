Project Juggler
===============

Concurrent event processor for NP experiments, based on the Gaudi framework.

Dependencies:
  - v1.8 requires Gaudi v33-34, ACTS 8.1+, DD4hep 1.16.1+ and eicd v0.2.0
  - v1.5 requires Gaudi v33-34, ACTS 8.1 and DD4hep 1.16.1

Overview
--------


### Internal Units

The juggler internal units are (`GeV`, `mm`, `ns`, and `radians`).

#### Units Table

| G4 | DD4hep | Gaudi | Juggler |
|----|--------|-------|---------|


### Running Juggler

Here is an example for topside
```
./scripts/run_topside.py -i inputevents.hepmc -o output_events.root -n 10
../where_ever/../juggler/build/run gaudirun.py options/example_reconstruction.py
```

## Components

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
