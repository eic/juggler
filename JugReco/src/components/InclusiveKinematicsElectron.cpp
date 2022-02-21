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

#include "eicd/VectorXYZT.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Reco {

class InclusiveKinematicsElectron : public GaudiAlgorithm {
public:
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticleCollection{
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

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton, m_neutron, m_electron;

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
    const eic::ReconstructedParticleCollection& parts = *(m_inputParticleCollection.get());
    // output collection
    auto& out_kinematics = *(m_outputInclusiveKinematicsCollection.createAndPut());

    // Loop over generated particles to get incoming electron beam
    eic::VectorXYZT ei, pi;
    bool found_electron = false;
    bool found_proton = false;
    //Also get the true scattered electron, which will not be included in the sum
    //over final-state particles for the JB reconstruction
    int32_t mcscatID = -1;
    for (const auto& p : mcparts) {
      if (p.getGeneratorStatus() == 4 && p.getPDG() == 11) {
        // Incoming electron
        ei.x = p.getMomentum().x;
        ei.y = p.getMomentum().y;
        ei.z = p.getMomentum().z;
        ei.t = p.getEnergy();

        //Should not include true event-by-event smearing of beam particles for reconstruction
        //Find a better way to do this...
        ei.x = 0.;
        ei.y = 0.;

        if( fabs(ei.z - 5.0) < 1.0 )
          ei.z = -5.0;
        else if( fabs(ei.z - 10.0) < 1.0 )
          ei.z = -10.0;
        else if( fabs(ei.z - 18.0) < 1.0 )
          ei.z = -18.0;
        
        ei.t = sqrt( ei.z*ei.z + m_electron*m_electron );

        found_electron = true;
      }
      if (p.getGeneratorStatus() == 4 && (p.getPDG() == 2212 || p.getPDG() == 2112)) {
        // Incoming proton
        pi.x = p.getMomentum().x;
        pi.y = p.getMomentum().y;
        pi.z = p.getMomentum().z;
        pi.t = p.getEnergy();

        //Should not include true event-by-event smearing of beam particles for reconstruction
        //Find a better way to do this...
        pi.y = 0;

        if( fabs(pi.z - 41.0) < 5.0 ){
          pi.x = 41.0*sin(m_crossingAngle);
          pi.z = 41.0*cos(m_crossingAngle);
        }
        else if( fabs(pi.z - 100.0) < 5.0 ){
          pi.x = 100.0*sin(m_crossingAngle);
          pi.z = 100.0*cos(m_crossingAngle);
        }
        else if( fabs(pi.z - 275.0) < 5.0 ){
          pi.x = 275.0*sin(m_crossingAngle);
          pi.z = 275.0*cos(m_crossingAngle);
        }
        pi.t = std::hypot(pi.x, pi.z, (p.getPDG() == 2212) ? m_proton : m_neutron);
      

        found_proton = true;
      }
      // Index of true Scattered electron. Currently taken as first status==1 electron in HEPMC record.
      if (p.getGeneratorStatus() == 1 && p.getPDG() == 11 && mcscatID == -1) {
        mcscatID = p.id();
      }

      if (found_electron && found_proton && mcscatID != -1) {
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
    if (mcscatID == -1) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No truth scattered electron found" << endmsg;
      }
      return StatusCode::SUCCESS;      
    }

    // Loop over reconstructed particles to get outgoing scattered electron
    // Use the true scattered electron from the MC information
    typedef std::pair<eic::VectorXYZT, eic::Index> t_electron;
    std::vector<t_electron> electrons;
    for (const auto& p : parts) {
      if (p.mcID().value == mcscatID) {
        // Outgoing electron
        electrons.push_back(t_electron(eic::VectorXYZT(p.p().x, p.p().y, p.p().z, p.energy()), eic::Index(p.ID())));
      }
    }

    if (electrons.size() > 0) {

      // DIS kinematics calculations
      auto kin = out_kinematics.create();
      const auto [ef, scatID] = electrons.front();
      const auto q = ei.subtract(ef);
      const auto q_dot_pi = q.dot(pi);
      kin.Q2(-q.dot(q));
      kin.y(q_dot_pi / ei.dot(pi));
      kin.nu(q_dot_pi / m_proton);
      kin.x(kin.Q2() / (2.*q_dot_pi));
      kin.W(sqrt(m_proton*m_proton + 2.*q_dot_pi - kin.Q2()));
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
