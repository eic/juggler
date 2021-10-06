#include <algorithm>
#include <cmath>
#include <fmt/format.h>

#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/VectorPolar.h"

namespace Jug::Fast {

class SmearedFarForwardParticles : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  DataHandle<dd4pod::Geant4ParticleCollection> m_inputParticles{"inputMCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"SmearedFarForwardParticles", Gaudi::DataHandle::Writer,
                                                                     this};

  Gaudi::Property<bool> m_enableZDC{this, "enableZDC", true};
  Gaudi::Property<bool> m_enableB0{this, "enableB0", true};
  Gaudi::Property<bool> m_enableRP{this, "enableRP", true};
  Gaudi::Property<bool> m_enableOMD{this, "enableOMD", true};

  // Beam energy, only used to determine the RP/OMD momentum ranges
  Gaudi::Property<double> m_ionBeamEnergy{this, "ionBeamEnergy", 100.};
  // RP default to 10-on-100 setting
  // Pz > 60% of beam energy (60% x 100GeV = 60GeV)
  // theta from 0.2mrad -> 5mrad
  Gaudi::Property<double> m_thetaMinRP{this, "thetaMinRP", 0.2e-3};
  Gaudi::Property<double> m_thetaMaxRP{this, "thetaMaxRP", 5e-3};
  Gaudi::Property<double> m_pMinRigidityRP{this, "pMinRigidityRP", 0.60};
  // B0
  Gaudi::Property<double> m_thetaMinB0{this, "thetaMinB0", 6.0e-3};
  Gaudi::Property<double> m_thetaMaxB0{this, "thetaMaxB0", 20.0e-3};
  // OMD default to 10-on-100 setting
  // 25% < P/Ebeam < 60% of beam energy (25% x 100GeV = 25GeV and 60% x 100GeV = 60GeV)
  // Angles both given for the small angle full-acceptance part,
  // and for the larger angle part where we only measure |phi| > rad
  Gaudi::Property<double> m_thetaMinFullOMD{this, "thetaMinFullOMD", 0.};
  Gaudi::Property<double> m_thetaMaxFullOMD{this, "thetaMaxFullOMD", 2e-3};
  Gaudi::Property<double> m_thetaMinPartialOMD{this, "thetaMinPartialOMD", 2.0e-3};
  Gaudi::Property<double> m_thetaMaxPartialOMD{this, "thetaMaxPartialOMD", 5.0e-3};
  Gaudi::Property<double> m_pMinRigidityOMD{this, "pMinRigidityOMD", 0.25};
  Gaudi::Property<double> m_pMaxRigidityOMD{this, "pMaxRigidityOMD", 0.60};

  // Crossing angle, set to -25mrad
  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025};

  Rndm::Numbers m_gaussDist;

  // Monte Carlo particle source identifier
  const int32_t m_kMonteCarloSource{uniqueID<int32_t>("mcparticles")};

  private:
  using RecData = eic::ReconstructedParticle;

  public:
  SmearedFarForwardParticles(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputParticles, "mcparticles");
    declareProperty("outputParticles", m_outputParticles, "ReconstructedParticles");
  }
  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;

    IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
    // use 0 for mean and 1 for standard deviation. Can rescale appropriately for the
    // different subsystems
    StatusCode sc = m_gaussDist.initialize(randSvc, Rndm::Gauss(0.0, 1.0));
    if (!sc.isSuccess()) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
  StatusCode execute() override {
    const auto& mc = *(m_inputParticles.get());
    auto& rc       = *(m_outputParticles.createAndPut());

    std::vector<std::vector<RecData>> rc_parts;
    if (m_enableZDC) {
      rc_parts.push_back(zdc(mc));
    }
    if (m_enableRP) {
      rc_parts.push_back(rp(mc));
    }
    if (m_enableB0) {
      rc_parts.push_back(b0(mc));
    }
    if (m_enableOMD) {
      rc_parts.push_back(omd(mc));
    }
    for (const auto& det : rc_parts) {
      for (const auto& part : det) {
        rc.push_back(part);
      }
    }
    return StatusCode::SUCCESS;
  }

private:
  // ZDC smearing as in eic_smear
  // https://github.com/eic/eicsmeardetectors/blob/9a1831dd97bf517b80a06043b9ee4bfb96b483d8/SmearMatrixDetector_0_1_FF.cxx#L224
  std::vector<RecData> zdc(const dd4pod::Geant4ParticleCollection& mc) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.genStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with genStatus = " << part.genStatus() << endmsg;
        }
        continue;
      }
      // only detect neutrons and photons
      const auto mom_ion = rotateLabToIonDirection(part.ps());
      if (part.pdgID() != 2112 && part.pdgID() != 22) {
        continue;
      }
      // only 0-->4.5 mrad
      if (mom_ion.theta() > 4.5 / 1000.) {
        continue;
      }
      const double E    = std::hypot(part.ps().mag(), part.mass());
      double conTerm = 0.05; //default 5%
      double stoTerm = 0.5;  //default 50%
      double angTerm = 0.003; //3mrad

      if(part.pdgID() == 2112){
        conTerm = 0.05; //default 5%
        stoTerm = 0.5;  //default 50%
        angTerm = 0.003; //3mrad
      }
      else if(part.pdgID() == 22){  //EMCAL expected to have slightly better performance
        conTerm = 0.03; //default 3%
        stoTerm = 0.25;  //default 25%
        angTerm = 0.003; //3mrad
      }

      const double dE   = sqrt((conTerm * E) * (conTerm * E) + stoTerm * stoTerm * E) * m_gaussDist(); //50%/SqrtE + 5%
      const double Es   = E + dE;
      const double th   = mom_ion.theta();
      const double dth  = (angTerm / sqrt(E)) * m_gaussDist();
      const double ths  = th + dth;
      const double phi  = mom_ion.phi();
      const double dphi = 0;
      const double phis = phi + dphi;
      const double moms = sqrt(Es * Es - part.mass() * part.mass());
      const eic::VectorPolar mom3s_ion{moms, ths, phis};
      const auto mom3s = rotateIonToLabDirection(mom3s_ion);
      eic::ReconstructedParticle rec_part;
      rec_part.ID({part.ID(), algorithmID()});
      rec_part.p(mom3s);
      rec_part.v({part.vs().x, part.vs().y, part.vs().z});
      rec_part.time(static_cast<float>(part.time()));
      rec_part.pid(part.pdgID());
      rec_part.status(0);
      rec_part.charge(static_cast<int16_t>(part.charge()));
      rec_part.weight(1.);
      rec_part.direction({mom3s.theta(), mom3s.phi()});
      rec_part.momentum(static_cast<float>(moms));
      rec_part.energy(static_cast<float>(Es));
      rec_part.mass(static_cast<float>(part.mass()));
      rec_part.mcID({part.ID(), m_kMonteCarloSource});
      rc.push_back(rec_part);

      if (msgLevel(MSG::DEBUG)) {
        debug()
            << fmt::format(
                   "Found ZDC particle: {}, Etrue: {}, Emeas: {}, ptrue: {}, pmeas: {}, theta_true: {}, theta_meas: {}",
                   part.pdgID(), E, rec_part.energy(), part.ps().mag(), rec_part.p().mag(), th, rec_part.p().theta())
            << endmsg;
      }
    }
    return rc;
  }
  // Fast B0 as in
  // https://github.com/eic/eicsmeardetectors/blob/9a1831dd97bf517b80a06043b9ee4bfb96b483d8/SmearMatrixDetector_0_1_FF.cxx#L254
  std::vector<RecData> b0(const dd4pod::Geant4ParticleCollection& mc) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.genStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with genStatus = " << part.genStatus() << endmsg;
        }
        continue;
      }
      // only detect charged hadrons and photons
      if (part.pdgID() != 2212 && part.pdgID() != -2212 && part.pdgID() != 211 && part.pdgID() != -211 &&
          part.pdgID() != 321 && part.pdgID() != -321 && part.pdgID() != 22) {
        continue;
      }
      // only 6-->20 mrad
      const auto mom_ion = rotateLabToIonDirection(part.ps());
      if (mom_ion.theta() < m_thetaMinB0 || mom_ion.theta() > m_thetaMaxB0) {
        continue;
      }
      rc.push_back(smearMomentum(part));
      if (msgLevel(MSG::DEBUG)) {
        auto& rec_part = rc.back();
        debug() << fmt::format("Found B0 particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.pdgID(), part.ps().mag(), rec_part.momentum(), std::hypot(part.ps().x, part.ps().y),
                               std::hypot(rec_part.p().x, rec_part.p().y), part.ps().theta(), rec_part.p().theta())
                << endmsg;
      }
    }

    return rc;
  }

  std::vector<RecData> rp(const dd4pod::Geant4ParticleCollection& mc) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.genStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with genStatus = " << part.genStatus() << endmsg;
        }
        continue;
      }
      // only detect protons
      if (part.pdgID() != 2212) {
        continue;
      }
      const auto mom_ion = rotateLabToIonDirection(part.ps());
      if (mom_ion.theta() < m_thetaMinRP || mom_ion.theta() > m_thetaMaxRP ||
          mom_ion.z < m_pMinRigidityRP * m_ionBeamEnergy) {
        continue;
      }
      rc.push_back(smearMomentum(part));
      if (msgLevel(MSG::DEBUG)) {
        auto& rec_part = rc.back();
        debug() << fmt::format("Found RP particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.pdgID(), part.ps().mag(), rec_part.momentum(), std::hypot(part.ps().x, part.ps().y),
                               std::hypot(rec_part.p().x, rec_part.p().y), part.ps().theta(), rec_part.p().theta())
                << endmsg;
      }
    }
    return rc;
  }

  std::vector<RecData> omd(const dd4pod::Geant4ParticleCollection& mc) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.genStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with genStatus = " << part.genStatus() << endmsg;
        }
        continue;
      }
      // only detect protons
      if (part.pdgID() != 2212) {
        continue;
      }
      const auto mom_ion = rotateLabToIonDirection(part.ps());
      if (mom_ion.z < m_pMinRigidityOMD * m_ionBeamEnergy || mom_ion.z > m_pMaxRigidityOMD * m_ionBeamEnergy) {
        continue;
      }
      // angle cut
      const double phi          = (mom_ion.phi() < M_PI) ? mom_ion.phi() : mom_ion.phi() - 2 * M_PI;
      const bool in_small_angle = (mom_ion.theta() > m_thetaMinFullOMD && mom_ion.theta() < m_thetaMaxFullOMD);
      const bool in_large_angle = (mom_ion.theta() > m_thetaMinPartialOMD && mom_ion.theta() < m_thetaMaxPartialOMD);
      if (!in_small_angle || (std::abs(phi) > 1 && !in_large_angle)) {
        continue;
      }
      rc.push_back(smearMomentum(part));
      if (msgLevel(MSG::DEBUG)) {
        auto& rec_part = rc.back();
        debug() << fmt::format("Found OMD particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.pdgID(), part.ps().mag(), rec_part.momentum(), std::hypot(part.ps().x, part.ps().y),
                               std::hypot(rec_part.p().x, rec_part.p().y), part.ps().theta(), rec_part.p().theta())
                << endmsg;
      }
    }
    return rc;
  }

  // all momentum smearing in EIC-smear for the far-forward region uses
  // the same 2 relations for P and Pt smearing (B0, RP, OMD)
  RecData smearMomentum(const dd4pod::ConstGeant4Particle& part) {
    const auto mom_ion = rotateLabToIonDirection(part.ps());
    const double p     = mom_ion.mag();
    const double dp    = (0.005 * p) * m_gaussDist();
    const double ps    = p + dp;

    const double pt  = std::hypot(mom_ion.x, mom_ion.y);
    const double dpt = (0.03 * pt) * m_gaussDist();
    // just apply relative smearing on px and py
    const double pxs = mom_ion.x + (1 + dpt / pt);
    const double pys = mom_ion.y + (1 + dpt / pt);
    // now get pz
    const double pzs = sqrt(ps * ps - pxs * pxs - pys * pys);
    // And build our 3-vector
    const eic::VectorXYZ psmear_ion = {pxs, pys, pzs};
    const auto psmear               = rotateIonToLabDirection(psmear_ion);
    eic::ReconstructedParticle rec_part;
    rec_part.ID({part.ID(), algorithmID()});
    rec_part.p(psmear);
    rec_part.v({part.vs().x, part.vs().y, part.vs().z});
    rec_part.time(static_cast<float>(part.time()));
    rec_part.pid(part.pdgID());
    rec_part.status(0);
    rec_part.charge(static_cast<int16_t>(part.charge()));
    rec_part.weight(1.);
    rec_part.direction({psmear.theta(), psmear.phi()});
    rec_part.momentum(static_cast<float>(ps));
    rec_part.energy(std::hypot(ps, part.mass()));
    rec_part.mass(static_cast<float>(part.mass()));
    rec_part.mcID({part.ID(), m_kMonteCarloSource});
    return rec_part;
  }

  // Rotate 25mrad about the y-axis
  eic::VectorXYZ rotateLabToIonDirection(const eic::VectorXYZ& vec) const {
    const double sth = sin(-m_crossingAngle);
    const double cth = cos(-m_crossingAngle);
    return {cth * vec.x + sth * vec.z, vec.y, -sth * vec.x + cth * vec.z};
  }
  eic::VectorXYZ rotateLabToIonDirection(const dd4pod::VectorXYZ& vec) const {
    return rotateLabToIonDirection(eic::VectorXYZ{vec.x, vec.y, vec.z});
  }

  eic::VectorXYZ rotateIonToLabDirection(const eic::VectorXYZ& vec) const {
    const double sth = sin(m_crossingAngle);
    const double cth = cos(m_crossingAngle);
    return {cth * vec.x + sth * vec.z, vec.y, -sth * vec.x + cth * vec.z};
  }
  eic::VectorXYZ rotateIonToLabDirection(const dd4pod::VectorXYZ& vec) const {
    return rotateIonToLabDirection(eic::VectorXYZ{vec.x, vec.y, vec.z});
  }

}; 

DECLARE_COMPONENT(SmearedFarForwardParticles)

} // namespace Jug::Fast

