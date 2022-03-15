#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/PhysicalConstants.h"
#include <algorithm>
#include <cmath>

#include "JugBase/IParticleSvc.h"
#include "JugBase/DataHandle.h"

#include "JugReco/Utilities/Beam.h"
#include "JugReco/Utilities/Boost.h"

#include "Math/GenVector/PxPyPzE4D.h"
typedef ROOT::Math::PxPyPzE4D<double> PxPyPzE4D;

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Reco {

class InclusiveKinematicsElectron : public GaudiAlgorithm {
private:
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticleCollection{
    "inputMCParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eicd::ReconstructedParticleCollection> m_inputParticleCollection{
    "ReconstructedParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eicd::InclusiveKinematicsCollection> m_outputInclusiveKinematicsCollection{
    "InclusiveKinematicsElectron",
    Gaudi::DataHandle::Writer,
    this};

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton, m_neutron, m_electron;

public:
  InclusiveKinematicsElectron(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "MCParticles");
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
    m_neutron = m_pidSvc->particle(2112).mass;
    m_electron = m_pidSvc->particle(11).mass;

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const edm4hep::MCParticleCollection& mcparts = *(m_inputMCParticleCollection.get());
    const eicd::ReconstructedParticleCollection& parts = *(m_inputParticleCollection.get());
    // output collection
    auto& out_kinematics = *(m_outputInclusiveKinematicsCollection.createAndPut());


    // Get incoming electron beam
    const auto ei_iter = Jug::Reco::Beam::find_first_beam_electron(mcparts);
    PxPyPzE4D ei;
    if (ei_iter != mcparts.end()) {
      ei = Jug::Reco::Beam::round_beam_four_momentum(
        ei_iter.getMomentum(),
        ei_iter.getEnergy(),
        m_electron,
        {-5.0, -10.0, -18.0},
        0.0);
    } else {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No beam electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    // Get incoming hadron beam
    const auto pi_iter = Jug::Reco::Beam::find_first_beam_hadron(mcparts);
    PxPyPzE4D pi;
    if (pi_iter != mcparts.end()) {
      pi = Jug::Reco::Beam::round_beam_four_momentum(
        pi_iter.getMomentum(),
        pi_iter.getEnergy(),
        pi_iter.getPDG() == 2212 ? m_proton : m_neutron,
        {41.0, 100.0, 275.0},
        m_crossingAngle);
    } else {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No beam hadron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    // Get first scattered electron
    const auto ef_iter = Jug::Reco::Beam::find_first_scattered_electron(mcparts);
    if (ef_iter == mcparts.end()) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No truth scattered electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    // Loop over reconstructed particles to get outgoing scattered electron
    // Use the true scattered electron from the MC information
    const auto rc_ef_iter = Jug::Reco::Beam::find_first_scattered_electron(parts);
    
    std::vector<PxPyPzE4D> electrons;
    auto is_electron = [](const auto &v) { return v.pdg == 11; };
    for (const auto& p: parts | std::views::filter(is_electron)) {
      electrons.push_back({p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy()});
    }
    std::range::sort(electrons, [](PxPyPzE4D a, PxPyPzE4D b) { return a.E() > b.E(); });

    if (electrons.size() > 0) {

      // DIS kinematics calculations
      auto kin = out_kinematics.create();
      const auto mass = pi_iter.getPDG() == 2212 ? m_proton : m_neutron;
      const auto ef = electrons.front();
      const auto q = ei.subtract(ef);
      const auto q_dot_pi = q.dot(pi);
      kin.setQ2(-q.dot(q));
      kin.setY(q_dot_pi / ei.dot(pi));
      kin.setNu(q_dot_pi / m_proton);
      kin.setX(kin.getQ2() / (2.*q_dot_pi));
      kin.setW(sqrt( + 2.*q_dot_pi - kin.Q2()));
      kin.setScat(part);

      // Debugging output
      if (msgLevel(MSG::DEBUG)) {
        debug() << "pi = " << pi << endmsg;
        debug() << "ei = " << ei << endmsg;
        debug() << "ef = " << ef << endmsg;
        debug() << "q = " << q << endmsg;
        debug() << "x,y,Q2,W,nu = "
                << kin.getX() << "," 
                << kin.getY() << ","
                << kin.getQ2() << ","
                << kin.getW() << ","
                << kin.getNu()
                << endmsg;
      }
    }

    return StatusCode::SUCCESS;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(InclusiveKinematicsElectron)

} // namespace Jug::Reco
