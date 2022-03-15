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

#include "JugReco/Utilities/Boost.h"

#include "Math/GenVector/PxPyPzE4D.h"
typedef ROOT::Math::PxPyPzE4D<double> PxPyPzE4D;

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Reco {

class InclusiveKinematicseSigma : public GaudiAlgorithm {
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
    "InclusiveKinematicseSigma",
    Gaudi::DataHandle::Writer,
    this};

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton, m_neutron, m_electron;

public:
  InclusiveKinematicseSigma(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "MCParticles");
    declareProperty("inputParticles", m_inputParticleCollection, "ReconstructedParticles");
    declareProperty("outputData", m_outputInclusiveKinematicsCollection, "InclusiveKinematicseSigma");
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
    eicd::VectorXYZT ei, pi;
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
        
        ei.t = sqrt( ei.z*ei.z + m_electron*m_electron);

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
    eicd::Index scatID;

    //Sums in colinear frame
    double pxsum = 0;
    double pysum = 0;
    double pzsum = 0;
    double Esum = 0;

    double pt_e = 0;
    double sigma_e = 0;

    for(const auto& p : parts){
        //Get the scattered electron index and angle
        if(p.mcID().value == mcscatID){

          TLorentzVector e_lab(p.p().x,p.p().y,p.p().z,p.energy());
          TLorentzVector e_boosted = Jug::Reco::Boost::apply_boost(ei,pi,e_lab);
          
          scatID = p.ID();
          pt_e = e_boosted.Pt();
          sigma_e = e_boosted.E() - e_boosted.Pz();
        }
        //Sum over all particles other than scattered electron
        else{
          
          //Lorentz vector in lab frame
          TLorentzVector hf_lab(p.p().x,p.p().y,p.p().z,p.energy());

          //Boost to colinear frame
          TLorentzVector hf_boosted = Jug::Reco::Boost::apply_boost(ei,pi,hf_lab);
          pxsum += hf_boosted.Px();
          pysum += hf_boosted.Py();
          pzsum += hf_boosted.Pz();
          Esum += hf_boosted.E();
        }
    }

    // DIS kinematics calculations
    auto sigma_h = Esum - pzsum;
    auto sigma_tot = sigma_e + sigma_h;

    if (scatID && sigma_h>0) {

      auto y_e = 1. - sigma_e / (2.*ei.energy());
      auto Q2_e = (pt_e*pt_e) / (1. - y_e);

      auto y_sig = sigma_h / sigma_tot;
      auto Q2_sig = (pt_e*pt_e) / (1. - y_sig);
      auto x_sig = Q2_sig / (4.*ei.energy()*pi.energy()*y_sig);
      
      auto Q2_esig = Q2_e;
      auto x_esig = x_sig;
      auto y_esig = Q2_esig / (4.*ei.energy()*pi.energy()*x_esig); //equivalent to (2*ei.energy() / sigma_tot)*y_sig 
      auto nu_esig = Q2_esig / (2.*m_proton*x_esig);
      auto W_esig = sqrt ( m_proton*m_proton + 2*m_proton*nu_esig - Q2_esig );      

      auto kin = out_kinematics.create();
      kin.setQ2(Q2_esig);
      kin.setY(y_esig);
      kin.setNu(nu_esig);
      kin.setX(x_esig);
      kin.setW(W_esig);
    
      kin.setScat(scatID); //MC index of scattered electron

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

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(InclusiveKinematicseSigma)

} // namespace Jug::Reco
