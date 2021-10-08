#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

#include <CLHEP/Vector/LorentzVector.h>

// Event Model related classes
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Fast {

class InclusiveKinematicsTruth : public GaudiAlgorithm, AlgorithmIDMixin<int32_t> {
public:
  DataHandle<dd4pod::Geant4ParticleCollection> m_inputParticleCollection{"mcparticles", Gaudi::DataHandle::Reader,
                                                                         this};
  DataHandle<eic::InclusiveKinematicsCollection> m_outputInclusiveKinematicsCollection{"InclusiveKinematicsTruth",
                                                                                       Gaudi::DataHandle::Writer, this};

  InclusiveKinematicsTruth(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputParticleCollection, "mcparticles");
    declareProperty("outputData", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsTruth");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collection
    const dd4pod::Geant4ParticleCollection& parts = *(m_inputParticleCollection.get());
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
    eic::Index scatID;

    for (const auto& p : parts) {

      if (p.genStatus() == 4 && p.pdgID() == 11) { // Incoming electron
        ei.setPx(p.ps().x);
        ei.setPy(p.ps().y);
        ei.setPz(p.ps().z);
        ei.setE(p.energy());
        ebeam_found = true;
      }
      if (p.genStatus() == 4 && p.pdgID() == 2122) { // Incoming proton
        pi.setPx(p.ps().x);
        pi.setPy(p.ps().y);
        pi.setPz(p.ps().z);
        pi.setE(p.energy());
        pbeam_found = true;
      }
      // Scattered electron. Currently taken as first status==1 electron in HEPMC record,
      // which seems to be correct based on a cursory glance at the Pythia8 output. In the future,
      // it may be better to trace back each final-state electron and see which one originates from
      // the beam.
      if (p.genStatus() == 1 && p.pdgID() == 11 && !scatID) {
        ef.setPx(p.ps().x);
        ef.setPy(p.ps().y);
        ef.setPz(p.ps().z);
        ef.setE(p.energy());

        scatID = p.ID();
      }
      if (ebeam_found && pbeam_found && scatID) {
        // all done!
        break;
      }
    }

    // DIS kinematics calculations
    auto kin = out_kinematics.create();
    const auto q    = ei - ef;
    kin.Q2(-1. * q.m2());
    kin.y((q * pi) / (ei * pi));
    kin.nu(q * pi / .938272);
    kin.x(kin.Q2() / (2. * q * pi));
    kin.W(sqrt((pi + q).m2()));
    kin.scatID(scatID);

    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(InclusiveKinematicsTruth)

} // namespace Jug::Fast
