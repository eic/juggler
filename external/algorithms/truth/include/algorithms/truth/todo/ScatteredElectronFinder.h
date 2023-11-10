// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/PhysicalConstants.h"
#include <algorithm>
#include <cmath>

#include "JugBase/IParticleSvc.h"
#include "JugBase/DataHandle.h"

#include "JugBase/Utilities/Beam.h"

#include "Math/Vector4D.h"
using ROOT::Math::PxPyPzEVector;

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"

namespace Jug::Fast {

class ScatteredElectronFinder : public GaudiAlgorithm {
private:
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticleCollection{
    "inputMCParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<edm4hep::MCParticleCollection> m_outputMCScatteredElectron{
    "outputMCScatteredElectron",
    Gaudi::DataHandle::Writer,
    this};

public:
  ScatteredElectronFinder(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "MCParticles");
    declareProperty("outputMCScatteredElectron", m_outputMCScatteredElectron, "MCScatteredElectron");
  }

  StatusCode initialize() override {
    return GaudiAlgorithm::initialize();
  }

  StatusCode execute() override {
    // input collections
    const auto& mcparts = *(m_inputMCParticleCollection.get());
    // output collection
    auto& out_electron = *(m_outputMCScatteredElectron.createAndPut());
    out_electron.setSubsetCollection();

    // Determine scattered electron
    //
    // Currently taken as first status==1 electron in HEPMC record,
    // which seems to be correct based on a cursory glance at the
    // Pythia8 output. In the future, it may be better to trace back
    // each final-state electron and see which one originates from
    // the beam.
    const auto ef_coll = Jug::Base::Beam::find_first_scattered_electron(mcparts);
    out_electron.push_back(ef_coll.front());

    return StatusCode::SUCCESS;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ScatteredElectronFinder)

} // namespace Jug::Fast
