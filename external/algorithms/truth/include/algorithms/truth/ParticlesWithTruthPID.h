// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

#include <algorithms/algorithm.h>

// Event Model related classes
#include <edm4eic/MCRecoParticleAssociationCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4hep/MCParticleCollection.h>

namespace algorithms::truth {

using ParticlesWithTruthPIDAlgorithm = Algorithm<
    Input<edm4hep::MCParticleCollection, edm4eic::TrackParametersCollection>,
    Output<edm4eic::ReconstructedParticleCollection, edm4eic::MCRecoParticleAssociationCollection>>;

class ParticlesWithTruthPID : public ParticlesWithTruthPIDAlgorithm {
public:
  ParticlesWithTruthPID(std::string_view name)
      : ParticlesWithTruthPIDAlgorithm{name,
                                       {"inputMCParticles", "inputTrackParameters"},
                                       {"outputParticles", "outputAssociations"},
                                       "Create mock reconstructed particles by associating truth "
                                       "PID with reconstructed tracks. Matching happens by "
                                       "comparing generated with reconstructed P, phi, and eta."} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  // Matching momentum tolerance requires 10% by default;
  Property<double> m_pRelativeTolerance{this, "pRelativeTolerance", 0.1,
                                        "Relative momentum tolarance factor"};
  // Matching phi tolerance of 30 mrad
  Property<double> m_phiTolerance{this, "phiTolerance", 0.030, "Azimuthal angle tolerance in rad"};
  // Matching eta tolerance of 0.2
  Property<double> m_etaTolerance{this, "etaTolerance", 0.2, "Eta tolerance"};
};

} // namespace algorithms::truth

