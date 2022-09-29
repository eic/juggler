// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Whitney Armstrong, Wouter Deconinck

#include <algorithms/algorithm.h>
#include <algorithms/random.h>

// Event Model related classes
#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4hep/MCParticleCollection.h>

namespace algorithms::truth {

using MC2SmearedParticleAlgorithm = Algorithm<Input<edm4hep::MCParticleCollection>,
                                              Output<edm4eic::ReconstructedParticleCollection>>;

class MC2SmearedParticle : public MC2SmearedParticleAlgorithm {
public:
  MC2SmearedParticle(std::string_view name)
      : MC2SmearedParticleAlgorithm{name,
                                    {"inputParticles"},
                                    {"outputParticles"},
                                    "Create mock reconstructed particles out of input MCParticles "
                                    "using Gaussian momentum smearing"} {}

  void init();
  void process(const Input&, const Output&) const;

private:
  Generator m_rng = RandomSvc::instance().generator();

  // (0.01 --> 1%)
  Property<double> m_smearing{this, "smearing", 0.01, "Sigma for Gaussian smearing factor"};
};

} // namespace algorithms::truth

