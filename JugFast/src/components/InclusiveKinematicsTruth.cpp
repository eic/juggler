#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/IParticleSvc.h"
#include "JugBase/DataHandle.h"

#include <CLHEP/Vector/LorentzVector.h>

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Fast {

class InclusiveKinematicsTruth : public GaudiAlgorithm {
public:
  DataHandle<edm4hep::MCParticleCollection> m_inputParticleCollection{"MCParticles", Gaudi::DataHandle::Reader,
                                                                         this};
  DataHandle<eicd::InclusiveKinematicsCollection> m_outputInclusiveKinematicsCollection{"InclusiveKinematicsTruth",
                                                                                       Gaudi::DataHandle::Writer, this};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton;
  double m_neutron;

  InclusiveKinematicsTruth(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputParticleCollection, "MCParticles");
    declareProperty("outputData", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsTruth");
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

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collection
    const edm4hep::MCParticleCollection& parts = *(m_inputParticleCollection.get());
    // output collection
    auto& out_kinematics = *(m_outputInclusiveKinematicsCollection.createAndPut());

    // Loop over generated particles to get incoming electron and proton beams
    // and the scattered electron. In the presence of QED radition on the incoming
    // or outgoing electron line, the vertex kinematics will be different than the
    // kinematics calculated using the scattered electron as done here.
    // Also need to update for CC events.
    auto ei = CLHEP::HepLorentzVector(0., 0., 0., 0.);
    auto pi = CLHEP::HepLorentzVector(0., 0., 0., 0.);
    auto ef = CLHEP::HepLorentzVector(0., 0., 0., 0.);

    bool ebeam_found = false;
    bool pbeam_found = false;
    int32_t mcscatID = -1;
  
    for (const auto& p : parts) {

      if (p.getGeneratorStatus() == 4 && p.getPDG() == 11) { // Incoming electron
        ei.setPx(p.getMomentum().x);
        ei.setPy(p.getMomentum().y);
        ei.setPz(p.getMomentum().z);
        ei.setE(p.getEnergy());
        ebeam_found = true;
      }
      else if (p.getGeneratorStatus() == 4 && p.getPDG() == 2212) { // Incoming proton
        pi.setPx(p.getMomentum().x);
        pi.setPy(p.getMomentum().y);
        pi.setPz(p.getMomentum().z);
        pi.setE(p.getEnergy());
        pbeam_found = true;
      }
      else if (p.getGeneratorStatus() == 4 && p.getPDG() == 2112) { // Incoming neutron
        pi.setPx(p.getMomentum().x);
        pi.setPy(p.getMomentum().y);
        pi.setPz(p.getMomentum().z);
        pi.setE(p.getEnergy());
        pbeam_found = true;
      }
      // Scattered electron. Currently taken as first status==1 electron in HEPMC record,
      // which seems to be correct based on a cursory glance at the Pythia8 output. In the future,
      // it may be better to trace back each final-state electron and see which one originates from
      // the beam.
      else if (p.getGeneratorStatus() == 1 && p.getPDG() == 11 && mcscatID == -1) {
        ef.setPx(p.getMomentum().x);
        ef.setPy(p.getMomentum().y);
        ef.setPz(p.getMomentum().z);
        ef.setE(p.getEnergy());

        mcscatID = p.id();
      }
      if (ebeam_found && pbeam_found && mcscatID != -1) {
        // all done!
        break;
      }
    }

    // Not all particles found
    if (ebeam_found == false) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No initial electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    if (pbeam_found == false) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No initial proton found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    if (mcscatID == -1) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No scattered electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    // DIS kinematics calculations
    auto kin = out_kinematics.create();
    const auto q = ei - ef;
    kin.setQ2(-1. * q.m2());
    kin.setY((q * pi) / (ei * pi));
    kin.setNu(q * pi / m_proton);
    kin.setX(kin.getQ2() / (2. * q * pi));
    kin.setW(sqrt((pi + q).m2()));
//    kin.scat(????); @TODO: this is now set as a OneToOneRelation to ReconstructedParticle, 
//                           which breaks for this algorithm that takes raw MCParticles

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

    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(InclusiveKinematicsTruth)

} // namespace Jug::Fast
