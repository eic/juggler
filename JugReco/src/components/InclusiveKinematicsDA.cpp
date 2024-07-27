// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Wouter Deconinck

#include <JugAlgo/Algorithm.h>
#include <EICrecon/algorithms/reco/InclusiveKinematicsDA.h>

namespace algorithms {

const std::shared_ptr<ParticleSvc::ParticleMap> ParticleSvc::kParticleMap =
  std::make_shared<ParticleSvc::ParticleMap>(ParticleSvc::ParticleMap{
    {           0, {           0,   0,   0.0           , "unknown" }},
  });

} // namespace algorithms

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
JUGALGO_DEFINE_ALGORITHM(InclusiveKinematicsDA, eicrecon::InclusiveKinematicsDA, Jug::Reco)
