// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

#include <algorithms/algorithm.h>
#include <algorithms/random.h>

// Event Model related classes
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/TrackParametersCollection.h"
#include "edm4hep/MCParticleCollection.h"

namespace algorithms::truth {

using ParticlesWithTruthPIDAlgorithm = Algorithm<
    Input<edm4hep::MCParticleCollection, edm4eic::TrackParametersCollection>,
    Output<edm4eic::ReconstructedParticleCollection, edm4eic::MCRecoParticleAssociationCollection>>;

class ParticlesWithTruthPID : public ParticlesWithTruthPIDAlgorithm {
public:
  ParticlesWithTruthPID(std::string_view name)
      : ParticlesWithTruthPIDAlgorithm{name,
                                      {"inputMCParticles", "inputTrackParameters"},
                                      {"outputParticles", "outputAssociations"}} {}

  void init();
  void process(const Input&, const Output&);

private:
  // Matching momentum tolerance requires 10% by default;
  Property<double> m_pRelativeTolerance{this, "pRelativeTolerance", {0.1}};
  // Matching phi tolerance of 10 mrad
  Property<double> m_phiTolerance{this, "phiTolerance", {0.030}};
  // Matchin eta tolerance of 0.1
  Property<double> m_etaTolerance{this, "etaTolerance", {0.2}};
};

} // namespace algorithms::truth

