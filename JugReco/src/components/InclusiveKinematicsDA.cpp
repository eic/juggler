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

class InclusiveKinematicsDA : public GaudiAlgorithm, AlgorithmIDMixin<int32_t> {
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
    "InclusiveKinematicsDA",
    Gaudi::DataHandle::Writer,
    this};

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton, m_electron;

  InclusiveKinematicsDA(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "mcparticles");
    declareProperty("inputParticles", m_inputParticleCollection, "ReconstructedParticles");
    declareProperty("outputData", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsDA");
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
    eic::Index mcscatID;
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
        
        ei.t = sqrt( ei.z*ei.z + m_electron*m_electron);

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
      if (p.genStatus() == 1 && p.pdgID() == 11 && !mcscatID) {
        mcscatID = p.ID();
      }

      if (found_electron && found_proton && mcscatID) {
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

    // Loop over reconstructed particles to get all outgoing particles
    // ----------------------------------------------------------------- 
    // Right now, everything is taken from Reconstructed particles branches.
    // 
    // This means the tracking detector is used for charged particles to caculate the momentum,
    // and the magnitude of this momentum plus the true PID to calculate the energy.
    // No requirement is made that these particles produce a hit in any other detector
    // 
    // Using the Reconstructed particles branches also means that the reconstruction for neutrals is done using the
    // calorimeter(s) information for the energy and angles, and then using this energy and the true PID to get the
    // magnitude of the momentum.
    // -----------------------------------------------------------------

    // Reconstructed Index of scattered electron
    eic::Index scatID;

    //Sums in colinear frame
    double pxsum = 0;
    double pysum = 0;
    double pzsum = 0;
    double Esum = 0;
    double theta_e = 0;
    for(const auto& p : parts){
        //Boost to colinear frame using first-order matrix
        double px_boosted = (-m_crossingAngle)/2.*p.energy() + p.p().x + (-m_crossingAngle)/2.*p.p().z;
        double py_boosted = p.p().y;
        double pz_boosted = (m_crossingAngle)/2.*p.p().x + p.p().z;
        double E_boosted = p.energy() + (-m_crossingAngle)/2.*p.p().x;

        //Get the scattered electron index and angle
        if(p.mcID() == mcscatID){
          scatID = p.ID();
          theta_e = acos( pz_boosted /sqrt(px_boosted*px_boosted + py_boosted*py_boosted + pz_boosted*pz_boosted) );
        }
        //Sum over all particles other than scattered electron
        else{
          pxsum += px_boosted;
          pysum += py_boosted;
          pzsum += pz_boosted;
          Esum += E_boosted;
        }
    }

    // DIS kinematics calculations
    auto sigma_h = Esum - pzsum;
    auto ptsum = sqrt(pxsum*pxsum + pysum*pysum);
    auto theta_h = 2.*atan(sigma_h/ptsum);

    if (scatID && sigma_h>0) {

      auto y_da = tan(theta_h/2.) / ( tan(theta_e/2.) + tan(theta_h/2.) );
      auto Q2_da = 4.*ei.energy()*ei.energy() * ( 1. / tan(theta_e/2.) ) * ( 1. / (tan(theta_e/2.) + tan(theta_h/2.)) );
      auto x_da = Q2_da / (4.*ei.energy()*pi.energy()*y_da);
      auto nu_da = Q2_da / (2.*m_proton*x_da);
      auto W_da = sqrt ( m_proton*m_proton + 2*m_proton*nu_da - Q2_da );      

      auto kin = out_kinematics.create();
      kin.Q2(Q2_da);
      kin.y(y_da);
      kin.nu(nu_da);
      kin.x(x_da);
      kin.W(W_da);
    
      kin.scatID(scatID); //MC index of scattered electron

      // Debugging output
      if (msgLevel(MSG::DEBUG)) {
        debug() << "pi = " << pi << endmsg;
        debug() << "ei = " << ei << endmsg;
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

DECLARE_COMPONENT(InclusiveKinematicsDA)

} // namespace Jug::Reco
