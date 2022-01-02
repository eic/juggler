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

#include "TLorentzVector.h"
#include "TVector3.h"

// Event Model related classes
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Reco {

class InclusiveKinematicsSigma : public GaudiAlgorithm, AlgorithmIDMixin<int32_t> {
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
    "InclusiveKinematicsSigma",
    Gaudi::DataHandle::Writer,
    this};

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton, m_neutron, m_electron;

  InclusiveKinematicsSigma(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "mcparticles");
    declareProperty("inputParticles", m_inputParticleCollection, "ReconstructedParticles");
    declareProperty("outputData", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsSigma");
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

  TLorentzVector apply_boost(eic::VectorXYZT ei_eicv, eic::VectorXYZT pi_eicv, TLorentzVector part){

    //Step 1: Find the needed boosts and rotations from the incoming lepton and hadron beams 
    //(note, this will give you a perfect boost, in principle you will not know the beam momenta exactly and should use an average)
  
    // Define the Boost to make beams back-to-back
    TLorentzVector ei(ei_eicv.x,ei_eicv.y,ei_eicv.z,ei_eicv.t);
    TLorentzVector pi(pi_eicv.x,pi_eicv.y,pi_eicv.z,pi_eicv.t);

    TLorentzVector cmBoost = (1./ei.E())*ei + (1./pi.E())*pi;

    TLorentzVector boost(-cmBoost.Px(),-cmBoost.Py(),-cmBoost.Pz(),cmBoost.E());
    TVector3 b;
    b = boost.BoostVector();

    TLorentzVector boostBack(0.0,0.0,cmBoost.Pz(),cmBoost.E());
    TVector3 bb;
    bb = boostBack.BoostVector(); // This will boost beams from a center of momentum frame back to (nearly) their original energies

    // Boost and rotate the incoming beams to find the proper rotations TLorentzVector
    pi.Boost(b); // Boost to COM frame
    ei.Boost(b);
    double rotAboutY = -1.0*TMath::ATan2(pi.Px(),pi.Pz()); // Rotate to remove x component of beams
    double rotAboutX = 1.0*TMath::ATan2(pi.Py(),pi.Pz()); // Rotate to remove y component of beams

    //Step 2: Apply boosts and rotations to any particle 4-vector 
    //(here too, choices will have to be made as to what the 4-vector is for reconstructed particles)
  
    //Boost and rotate particle 4-momenta into the headon frame
    part.Boost(b);
    part.RotateY(rotAboutY);
    part.RotateX(rotAboutX);
    part.Boost(bb);

    return part;

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
        
        ei.t = sqrt( ei.z*ei.z + m_electron*m_electron);

        found_electron = true;
      }
      if (p.genStatus() == 4 && (p.pdgID() == 2212 || p.pdgID() == 2112)) {
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

        pi.t = std::hypot(pi.x, pi.z, (p.pdgID() == 2212) ? m_proton : m_neutron);

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

    double pt_e = 0;
    double sigma_e = 0;

    for(const auto& p : parts){
        //Get the scattered electron index and angle
        if(p.mcID().value == mcscatID){

          TLorentzVector e_lab(p.p().x,p.p().y,p.p().z,p.energy());
          TLorentzVector e_boosted = apply_boost(ei,pi,e_lab);
          
          scatID = p.ID();
          pt_e = e_boosted.Pt();
          sigma_e = e_boosted.E() - e_boosted.Pz();
        }
        //Sum over all particles other than scattered electron
        else{
          
          //Lorentz vector in lab frame
          TLorentzVector hf_lab(p.p().x,p.p().y,p.p().z,p.energy());

          //Boost to colinear frame
          TLorentzVector hf_boosted = apply_boost(ei,pi,hf_lab);
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

      auto y_sig = sigma_h / sigma_tot;
      auto Q2_sig = (pt_e*pt_e) / (1. - y_sig);
      auto x_sig = Q2_sig / (4.*ei.energy()*pi.energy()*y_sig);
      auto nu_sig = Q2_sig / (2.*m_proton*x_sig);
      auto W_sig = sqrt ( m_proton*m_proton + 2*m_proton*nu_sig - Q2_sig );      

      auto kin = out_kinematics.create();
      kin.Q2(Q2_sig);
      kin.y(y_sig);
      kin.nu(nu_sig);
      kin.x(x_sig);
      kin.W(W_sig);
    
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

DECLARE_COMPONENT(InclusiveKinematicsSigma)

} // namespace Jug::Reco
