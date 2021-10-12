#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

#include "eicd/VectorXYZT.h"

// Event Model related classes
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Reco {

class InclusiveKinematicsElectron : public GaudiAlgorithm, AlgorithmIDMixin<int32_t> {
public:
  DataHandle<dd4pod::Geant4ParticleCollection> m_inputMCParticleCollection{
    "inputMCParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eic::ReconstructedParticleCollection> m_inputParticleCollection{
    "ReconstructedParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eic::InclusiveKinematicsCollection> m_outputInclusiveKinematicsCollection{
    "InclusiveKinematicsElectron",
    Gaudi::DataHandle::Writer,
    this};

  InclusiveKinematicsElectron(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "mcparticles");
    declareProperty("inputParticles", m_inputParticleCollection, "ReconstructedParticles");
    declareProperty("outputData", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsElectron");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const dd4pod::Geant4ParticleCollection& mcparts = *(m_inputMCParticleCollection.get());
    const eic::ReconstructedParticleCollection& parts = *(m_inputParticleCollection.get());
    // output collection
    auto& out_kinematics = *(m_outputInclusiveKinematicsCollection.createAndPut());

    // Loop over generated particles to get incoming electron beam
    eic::VectorXYZT ei, pi;
    bool found_electron = false;
    bool found_proton = false;
    for (const auto& p : mcparts) {
      if (p.genStatus() == 4 && p.pdgID() == 11) {
        // Incoming electron
        ei.x = p.ps().x;
        ei.y = p.ps().y;
        ei.z = p.ps().z;
        ei.t = p.energy();
        found_electron = true;
      }
      if (p.genStatus() == 4 && p.pdgID() == 2212) {
        // Incoming proton
        pi.x = p.ps().x;
        pi.y = p.ps().y;
        pi.z = p.ps().z;
        pi.t = p.energy();
        found_proton = true;
      }
      if (found_electron && found_proton) {
        break;
      }
    }

    // Loop over reconstructed particles to get outgoing electrons
    std::vector<std::pair<eic::VectorXYZT, eic::Index>> ef;
    for (const auto& p : parts) {
      if (p.pid() == 11) {
        // Outgoing electron
        ef.push_back(std::make_pair<eic::VectorXYZT, eic::Index>(eic::VectorXYZT(p.p().x, p.p().y, p.p().z, p.energy()), eic::Index(p.ID())));
      }
    }

    if (ef.size() > 0) {
      // Sort by momentum magnitude
      std::sort(ef.begin(), ef.end(), [](std::pair<eic::VectorXYZT, eic::Index> a,  std::pair<eic::VectorXYZT, eic::Index> b) {
        return a.first.mag() > b.first.mag();
      });

      // DIS kinematics calculations
      auto kin = out_kinematics.create();
      const auto q = ei.subtract(ef.front().first);
      kin.Q2(q.mass());
      kin.y(q.dot(pi) / ei.dot(pi));
      kin.nu(q.dot(pi) / .938272);
      kin.x(kin.Q2() / (2. * q.dot(pi)));
      kin.W(sqrt(pi.add(q).mass()));
      kin.scatID(ef.front().second);
    }

    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(InclusiveKinematicsElectron)

} // namespace Jug::Reco
