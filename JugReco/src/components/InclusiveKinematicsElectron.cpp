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

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton, m_electron;

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
    m_electron = m_pidSvc->particle(11).mass;

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
    //Also get the true scattered electron, which will not be included in the sum
    //over final-state particles for the JB reconstruction
    int32_t mcscatID = -1;
    for (const auto& p : mcparts) {
      if (p.genStatus() == 4 && p.pdgID() == 11) {
        // Incoming electron
        ei.x = p.ps().x;
        ei.y = p.ps().y;
        ei.z = p.ps().z;
        ei.t = p.energy();

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
      if (p.genStatus() == 4 && p.pdgID() == 2212) {
        // Incoming proton
        pi.x = p.ps().x;
        pi.y = p.ps().y;
        pi.z = p.ps().z;
        pi.t = p.energy();

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

        pi.t = sqrt( pi.x*pi.x + pi.z*pi.z + m_proton*m_proton);

        found_proton = true;
      }
      // Index of true Scattered electron. Currently taken as first status==1 electron in HEPMC record.
      if (p.genStatus() == 1 && p.pdgID() == 11 && mcscatID == -1) {
        mcscatID = p.ID();
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
