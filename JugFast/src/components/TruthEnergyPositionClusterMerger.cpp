// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Wouter Deconinck

#include <JugAlgo/Algorithm.h>
#include <DD4hep/DD4hepUnits.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <EICrecon/algorithms/calorimetry/TruthEnergyPositionClusterMerger.h>

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
JUGALGO_DEFINE_ALGORITHM(TruthEnergyPositionClusterMerger, eicrecon::TruthEnergyPositionClusterMerger, Jug::Fast)
