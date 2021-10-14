#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/IParticleSvc.h"
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

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton;

  InclusiveKinematicsElectron(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "mcparticles");
    declareProperty("inputParticles", m_inputParticleCollection, "ReconstructedParticles");
    declareProperty("outputData", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsElectron");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;

    m_pidSvc = service("ParticleSvc");
    if (!m_pidSvc) {
      error() << "Unable to locate Particle Service. "
              << "Make sure you have ParticleSvc in the configuration."
              << endmsg;
      return StatusCode::FAILURE;
    }
    m_proton = m_pidSvc->particle(2212).mass;

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
    if (found_electron == false) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No initial electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    if (found_proton == false) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No initial proton found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    // Loop over reconstructed particles to get outgoing electrons
    typedef std::pair<eic::VectorXYZT, eic::Index> t_electron;
    std::vector<t_electron> electrons;
    for (const auto& p : parts) {
      if (p.pid() == 11) {
        // Outgoing electron
        electrons.push_back(t_electron(eic::VectorXYZT(p.p().x, p.p().y, p.p().z, p.energy()), eic::Index(p.ID())));
      }
    }

    if (electrons.size() > 0) {
      // Sort by momentum magnitude
      auto compare = [](const t_electron& a, const t_electron& b) {
        return a.first.mag() > b.first.mag();
      };
      std::sort(electrons.begin(), electrons.end(), compare);

      // DIS kinematics calculations
      auto kin = out_kinematics.create();
      const auto [ef, scatID] = electrons.front();
      const auto q = ei.subtract(ef);
      const auto q_dot_pi = q.dot(pi);
      kin.Q2(-q.dot(q));
      kin.y(q_dot_pi / ei.dot(pi));
      kin.nu(q_dot_pi / m_proton);
      kin.x(kin.Q2() / (2.*q_dot_pi));
      kin.W(m_proton*m_proton + 2.*q_dot_pi - kin.Q2());
      kin.scatID(scatID);

      // Debugging output
      if (msgLevel(MSG::DEBUG)) {
        debug() << "pi = " << pi << endmsg;
        debug() << "ei = " << ei << endmsg;
        debug() << "ef = " << ef << endmsg;
        debug() << "q = " << q << endmsg;
        debug() << "x,y,Q2,W,nu = "
                << kin.x() << "," 
                << kin.y() << ","
                << kin.Q2() << ","
                << kin.W() << ","
                << kin.nu()
                << endmsg;
      }
    }

    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(InclusiveKinematicsElectron)

} // namespace Jug::Reco
