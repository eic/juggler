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

    // Loop over generated particles to get incoming electron beam
    PxPyPzE4D ei, pi;
    bool found_electron = false;
    bool found_proton = false;
    //Also get the true scattered electron, which will not be included in the sum
    //over final-state particles for the JB reconstruction
    auto ef_iter = mcparts.end();
    for (const auto& p : mcparts) {
      if (p.getGeneratorStatus() == 4 && p.getPDG() == 11) {
        // Incoming electron
        ei.SetPxPyPzE(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());

        //Should not include true event-by-event smearing of beam particles for reconstruction
        //Find a better way to do this...
        ei.SetPx(0.);
        ei.SetPy(0.);

        if( fabs(ei.Pz() - 5.0) < 1.0 )
          ei.SetPz(-5.0);
        else if( fabs(ei.Pz() - 10.0) < 1.0 )
          ei.SetPz(-10.0);
        else if( fabs(ei.Pz() - 18.0) < 1.0 )
          ei.SetPz(-18.0);
        
        ei.SetE(std::hypot(ei.Pz(), m_electron));

        found_electron = true;
      }
      if (p.getGeneratorStatus() == 4 && (p.getPDG() == 2212 || p.getPDG() == 2112)) {
        // Incoming proton
        pi.SetPxPyPzE(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());

        //Should not include true event-by-event smearing of beam particles for reconstruction
        //Find a better way to do this...
        pi.SetPy(0.);

        if( fabs(pi.Pz() - 41.0) < 5.0 ){
          pi.SetPx(41.0*sin(m_crossingAngle));
          pi.SetPz(41.0*cos(m_crossingAngle));
        }
        else if( fabs(pi.Pz() - 100.0) < 5.0 ){
          pi.SetPx(100.0*sin(m_crossingAngle));
          pi.SetPz(100.0*cos(m_crossingAngle));
        }
        else if( fabs(pi.Pz() - 275.0) < 5.0 ){
          pi.SetPx(275.0*sin(m_crossingAngle));
          pi.SetPz(275.0*cos(m_crossingAngle));
        }
        pi.SetE(std::hypot(pi.Px(), pi.Pz(), (p.getPDG() == 2212) ? m_proton : m_neutron));
      
        found_proton = true;
      }
      // Index of true Scattered electron. Currently taken as first status==1 electron in HEPMC record.
      if (p.getGeneratorStatus() == 1 && p.getPDG() == 11 && ef_iter == mcparts.end()) {
        ef_iter = p;
      }

      if (found_electron && found_proton && ef_iter != mcparts.end()) {
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
    if (ep == mcparts.end()) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No truth scattered electron found" << endmsg;
      }
      return StatusCode::SUCCESS;      
    }

    // Loop over reconstructed particles to get outgoing scattered electron
    // Use the true scattered electron from the MC information
    typedef std::pair<PxPyPzE4D, eic::ReconstructedParticleCollection::iterator> t_electron;
    std::vector<t_electron> electrons;
    for (const auto& p: parts) {
      if (p.mcID() == ef_iter) {
        // Outgoing electron
        electrons.push_back(t_electron(PxPyPzE4D(p.p().x, p.p().y, p.p().z, p.energy()), p));
      }
    }

    if (electrons.size() > 0) {

      // DIS kinematics calculations
      auto kin = out_kinematics.create();
      const auto [ef, part] = electrons.front();
      const auto q = ei.subtract(ef);
      const auto q_dot_pi = q.dot(pi);
      kin.Q2(-q.dot(q));
      kin.y(q_dot_pi / ei.dot(pi));
      kin.nu(q_dot_pi / m_proton);
      kin.x(kin.Q2() / (2.*q_dot_pi));
      kin.W(sqrt(m_proton*m_proton + 2.*q_dot_pi - kin.Q2()));
      kin.scatID(part);

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

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(InclusiveKinematicsElectron)

} // namespace Jug::Reco
