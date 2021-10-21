[TOC]

# Project Juggler

Concurrent event processor for NP experiments, based on the Gaudi framework.

Dependencies:
  - v4.x requires Gaudi v36+, ACTS v13+, DD4hep 1.17+, NPdet v1.0+ and eicd v1.1+
  - v3.6 requires Gaudi v36+, ACTS v13+, DD4hep 1.17+, NPdet v1.0.0 and eicd v0.9.0
  - v3.5 requires Gaudi v36+, ACTS v13+, DD4hep 1.17+, NPdet v0.9.0, eicd v0.8.0
  - v3.4 requires Gaudi v36+, ACTS v8.2+, DD4hep 1.17+, NPdet v0.9.0, eicd v0.8.0
  - v3.0 requires Gaudi v35+, ACTS v8.1+, DD4hep 1.17+, NPdet v0.7.0 eicd v0.5.0
  - v2.0 requires Gaudi v35+, ACTS v8.1+, DD4hep 1.17+ and eicd v0.2.0
  - v1.8 requires Gaudi v33-34, ACTS 8.1+, DD4hep 1.16.1+ and eicd v0.2.0
  - v1.5 requires Gaudi v33-34, ACTS 8.1 and DD4hep 1.16.1

Overview
--------

### Components

 - [Juggler Base](@ref base)
 - [Digitization Algorithms](@ref digi)
 - [Hit Reconstruction Algorithms](@ref reco)
 - [Tracking Algorithms](@ref tracking) 
 - [PID Algorithms](@ref pid) 

### Internal Units

The juggler internal units are (`GeV`, `mm`, `ns`, and `radians`).


### Running Juggler

Here is an example for topside
```
./scripts/run_topside.py -i inputevents.hepmc -o output_events.root -n 10
../where_ever/../juggler/build/run gaudirun.py options/example_reconstruction.py
```

# Outline of tracking and vertexing

## The ACTS way of tracking

First, the geometry has to be constructed.  Assuming this is already done, we describe the data processing.

### Source Links and Measurements

```
  Tracker Hit  +  Geometry --> Measurement
  Hit CellID   +  Surface  --> Source Link
```

A `SourceLinker` is an algorithm that looks up the surface and maps it to the hit's `cellID`.
Naturally the same algorithm outputs measurements which are also mapped to the hit and contain the position 
and sensor size information (via covariance matrix).


### Proto tracks

Both track finding and fitting will use the information contained in the source links and measurements.
Track finding produces `proto tracks` or groupings of hits.  Each proto track is simply a `std::vector<int>` storing the index of the
hits associated with a track seed.

### Initial Track parameters and Seeding

A Kalman filter needs a starting point and those are the `Initial Track Parameters`. These can be determined many different ways. 
Conceptually the process of determining these parameters begins with track seeding.


## Vertexing

WIP





