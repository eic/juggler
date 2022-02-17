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
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/VectorPolar.h"

namespace {
  enum DetectorTags {
    kTagB0 = 1,
    kTagRP = 2,
    kTagOMD = 3,
    kTagZDC = 4
  };
}

namespace Jug::Fast {

class SmearedFarForwardParticles : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"inputMCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"SmearedFarForwardParticles", Gaudi::DataHandle::Writer,
                                                                     this};

  Gaudi::Property<bool> m_enableZDC{this, "enableZDC", true};
  Gaudi::Property<bool> m_enableB0{this, "enableB0", true};
  Gaudi::Property<bool> m_enableRP{this, "enableRP", true};
  Gaudi::Property<bool> m_enableOMD{this, "enableOMD", true};

  // Beam energy, only used to determine the RP/OMD momentum ranges
  Gaudi::Property<float> m_ionBeamEnergy{this, "ionBeamEnergy", 0.};
  // RP default to 10-on-100 setting
  // Pz > 60% of beam energy (60% x 100GeV = 60GeV)
  // theta from 0.2mrad -> 5mrad
  Gaudi::Property<float> m_thetaMinRP{this, "thetaMinRP", 0.2e-3};
  Gaudi::Property<float> m_thetaMaxRP{this, "thetaMaxRP", 5e-3};
  Gaudi::Property<float> m_pMinRigidityRP{this, "pMinRigidityRP", 0.60};
  // B0
  Gaudi::Property<float> m_thetaMinB0{this, "thetaMinB0", 6.0e-3};
  Gaudi::Property<float> m_thetaMaxB0{this, "thetaMaxB0", 20.0e-3};
  // OMD default to 10-on-100 setting
  // 25% < P/Ebeam < 60% of beam energy (25% x 100GeV = 25GeV and 60% x 100GeV = 60GeV)
  // Angles both given for the small angle full-acceptance part,
  // and for the larger angle part where we only measure |phi| > rad
  Gaudi::Property<float> m_thetaMinFullOMD{this, "thetaMinFullOMD", 0.};
  Gaudi::Property<float> m_thetaMaxFullOMD{this, "thetaMaxFullOMD", 2e-3};
  Gaudi::Property<float> m_thetaMinPartialOMD{this, "thetaMinPartialOMD", 2.0e-3};
  Gaudi::Property<float> m_thetaMaxPartialOMD{this, "thetaMaxPartialOMD", 5.0e-3};
  Gaudi::Property<float> m_pMinRigidityOMD{this, "pMinRigidityOMD", 0.25};
  Gaudi::Property<float> m_pMaxRigidityOMD{this, "pMaxRigidityOMD", 0.60};

  // Crossing angle, set to -25mrad
  Gaudi::Property<float> m_crossingAngle{this, "crossingAngle", -0.025}; //-0.025}; -- causes double rotation with afterburner

  Rndm::Numbers m_gaussDist;

  // Monte Carlo particle source identifier
  const int32_t m_kMonteCarloSource{uniqueID<int32_t>("MCParticles")};

  private:
  using RecData = eic::ReconstructedParticle;

  public:
  SmearedFarForwardParticles(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputMCParticles, "MCParticles");
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
    const auto& mc = *(m_inputMCParticles.get());
    auto& rc       = *(m_outputParticles.createAndPut());

    double ionBeamEnergy = 0;
    if (m_ionBeamEnergy > 0) {
      ionBeamEnergy = m_ionBeamEnergy;
    } else {
      for (const auto& part: mc) {
        if (part.getGeneratorStatus() == 4 && part.getPDG() == 2212) {
          auto E = part.getEnergy();
          if (33 < E && E < 50) {
            ionBeamEnergy = 41;
          } else if (80 < E && E < 120) {
            ionBeamEnergy = 100;
          } else if (220 < E && E < 330) {
            ionBeamEnergy = 275;
          } else {
            warning() << "Ion beam energy " << E << " not a standard setting." << endmsg;
            ionBeamEnergy = E;
          }
          break;
        }
      }
      if (ionBeamEnergy == 0) {
        warning() << "No incoming ion beam; using 100 GeV ion beam energy." << endmsg;
        ionBeamEnergy = 100;
      }
    }

    std::vector<std::vector<RecData>> rc_parts;
    if (m_enableZDC) {
      rc_parts.push_back(zdc(mc, ionBeamEnergy));
    }
    if (m_enableRP) {
      rc_parts.push_back(rp(mc, ionBeamEnergy));
    }
    if (m_enableB0) {
      rc_parts.push_back(b0(mc, ionBeamEnergy));
    }
    if (m_enableOMD) {
      rc_parts.push_back(omd(mc, ionBeamEnergy));
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
  std::vector<RecData> zdc(const edm4hep::MCParticleCollection& mc, const double ionBeamEnergy) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.getGeneratorStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with getGeneratorStatus = " << part.getGeneratorStatus() << endmsg;
        }
        continue;
      }
      // only detect neutrons and photons
      const auto mom_ion = rotateLabToIonDirection(part.getMomentum());
      if (part.getPDG() != 2112 && part.getPDG() != 22) {
        continue;
      }
      // only 0-->4.5 mrad
      const auto mom_ion_phi = std::atan2(mom_ion.x, mom_ion.y);
      const auto mom_ion_theta = std::atan2(std::hypot(mom_ion.x, mom_ion.y), mom_ion.z);
      if (mom_ion_theta > 4.5 / 1000.) {
        continue;
      }
      const auto E = part.getEnergy();
      float conTerm = 0.05; //default 5%
      float stoTerm = 0.5;  //default 50%
      float angTerm = 0.003; //3mrad

      if(part.getPDG() == 2112){
        conTerm = 0.05; //default 5%
        stoTerm = 0.5;  //default 50%
        angTerm = 0.003; //3mrad
      }
      else if(part.getPDG() == 22){  //EMCAL expected to have slightly better performance
        conTerm = 0.03; //default 3%
        stoTerm = 0.10;  //default 10% for WSciFi
        angTerm = 0.001; //1mrad is the detault for the block size
      }

      const auto dE   = sqrt((conTerm * E) * (conTerm * E) + stoTerm * stoTerm * E) * m_gaussDist(); //50%/SqrtE + 5%
      const auto Es   = E + dE;
      const auto th   = mom_ion_theta;
      const auto dth  = (angTerm / sqrt(E)) * m_gaussDist();
      const auto ths  = th + dth;
      const auto phi  = mom_ion_phi;
      const auto dphi = 0;
      const auto phis = phi + dphi;
      const auto moms = sqrt(Es * Es - part.getMass() * part.getMass());
      const edm4hep::Vector3f mom3s_ion{moms*sin(ths)*cos(phis), moms*sin(ths)*sin(phis), moms*cos(ths)};
      const auto mom3s = rotateIonToLabDirection(mom3s_ion);
      const auto mom3s_phi = std::atan2(mom3s.y, mom3s.x);
      const auto mom3s_theta = std::atan2(std::hypot(mom3s.y, mom3s.x), mom3s.z);
      eic::ReconstructedParticle rec_part;
      rec_part.ID({part.id(), algorithmID()});
      rec_part.p({mom3s.x, mom3s.y, mom3s.z});
      rec_part.v({part.getVertex().x, part.getVertex().y, part.getVertex().z});
      rec_part.time(static_cast<float>(part.getTime()));
      rec_part.pid(part.getPDG());
      rec_part.status(kTagZDC);
      rec_part.charge(static_cast<int16_t>(part.getCharge()));
      rec_part.weight(1.);
      rec_part.direction({mom3s_theta, mom3s_phi});
      rec_part.momentum(static_cast<float>(moms));
      rec_part.energy(static_cast<float>(Es));
      rec_part.mass(static_cast<float>(part.getMass()));
      rc.push_back(rec_part);

      if (msgLevel(MSG::DEBUG)) {
        const auto& part_p = part.getMomentum();
        const auto part_p_mag = std::hypot(part_p.x, part_p.y, part_p.z);
        debug()
            << fmt::format(
                   "Found ZDC particle: {}, Etrue: {}, Emeas: {}, ptrue: {}, pmeas: {}, theta_true: {}, theta_meas: {}",
                   part.getPDG(), E, rec_part.energy(), part_p_mag, rec_part.p().mag(), th, rec_part.p().theta())
            << endmsg;
      }
    }
    return rc;
  }
  // Fast B0 as in
  // https://github.com/eic/eicsmeardetectors/blob/9a1831dd97bf517b80a06043b9ee4bfb96b483d8/SmearMatrixDetector_0_1_FF.cxx#L254
  std::vector<RecData> b0(const edm4hep::MCParticleCollection& mc, const double ionBeamEnergy) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.getGeneratorStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with getGeneratorStatus = " << part.getGeneratorStatus() << endmsg;
        }
        continue;
      }
      // only detect charged hadrons and photons
      if (part.getPDG() != 2212 && part.getPDG() != -2212 && part.getPDG() != 211 && part.getPDG() != -211 &&
          part.getPDG() != 321 && part.getPDG() != -321 && part.getPDG() != 22) {
        continue;
      }
      // only 6-->20 mrad
      const auto mom_ion = removeCrossingAngle(part.getMomentum()); //rotateLabToIonDirection(part.getMomentum());
      const auto mom_ion_theta = std::atan2(std::hypot(mom_ion.x, mom_ion.y), mom_ion.z);
      if (mom_ion_theta < m_thetaMinB0 || mom_ion_theta > m_thetaMaxB0) {
        continue;
      }
      auto rc_part = smearMomentum(part);
      // we don't detect photon energy, just its angles and presence
      if (part.getPDG() == 22) {
        rc_part.p({0,0,0});
        rc_part.energy(0);
      }
      rc_part.status(kTagB0);
      rc.push_back(rc_part);
      if (msgLevel(MSG::DEBUG)) {
        auto& rec_part = rc.back();
        const auto& part_p = part.getMomentum();
        const auto part_p_pt = std::hypot(part_p.x, part_p.y);
        const auto part_p_mag = std::hypot(part_p.x, part_p.y, part_p.z);
        const auto part_p_theta = std::atan2(std::hypot(part_p.x, part_p.y), part_p.z);
        debug() << fmt::format("Found B0 particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.getPDG(), part_p_mag, rec_part.momentum(), part_p_pt,
                               std::hypot(rec_part.p().x, rec_part.p().y), part_p_theta, rec_part.p().theta())
                << endmsg;
      }
    }

    return rc;
  }

  std::vector<RecData> rp(const edm4hep::MCParticleCollection& mc, const double ionBeamEnergy) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.getGeneratorStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with getGeneratorStatus = " << part.getGeneratorStatus() << endmsg;
        }
        continue;
      }
      // only detect protons
      if (part.getPDG() != 2212) {
        continue;
      }
      const auto mom_ion = removeCrossingAngle(part.getMomentum()); //rotateLabToIonDirection(part.getMomentum());
      const auto mom_ion_theta = std::atan2(std::hypot(mom_ion.x, mom_ion.y), mom_ion.z);
      if (mom_ion_theta < m_thetaMinRP || mom_ion_theta > m_thetaMaxRP ||
          mom_ion.z < m_pMinRigidityRP * ionBeamEnergy) {
        continue;
      }
      auto rc_part = smearMomentum(part);
      rc_part.status(kTagRP);
      rc.push_back(rc_part);
      if (msgLevel(MSG::DEBUG)) {
        auto& rec_part = rc.back();
        const auto& part_p = part.getMomentum();
        const auto part_p_pt = std::hypot(part_p.x, part_p.y);
        const auto part_p_mag = std::hypot(part_p.x, part_p.y, part_p.z);
        const auto part_p_theta = std::atan2(std::hypot(part_p.x, part_p.y), part_p.z);
        debug() << fmt::format("Found RP particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.getPDG(), part_p_mag, rec_part.momentum(), part_p_pt,
                               std::hypot(rec_part.p().x, rec_part.p().y), part_p_theta, rec_part.p().theta())
                << endmsg;
      }
    }
    return rc;
  }

  std::vector<RecData> omd(const edm4hep::MCParticleCollection& mc, const double ionBeamEnergy) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.getGeneratorStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with getGeneratorStatus = " << part.getGeneratorStatus() << endmsg;
        }
        continue;
      }
      // only detect protons
      if (part.getPDG() != 2212) {
        continue;
      }
      const auto mom_ion = removeCrossingAngle(part.getMomentum()); //rotateLabToIonDirection(part.getMomentum());
      //const auto mom_ion_phi = std::atan2(mom_ion.x, mom_ion.y);
      //const auto mom_ion_theta = std::atan2(std::hypot(mom_ion.x, mom_ion.y), mom_ion.z);
      if (mom_ion.z < m_pMinRigidityOMD * ionBeamEnergy || mom_ion.z > m_pMaxRigidityOMD * ionBeamEnergy) {
        continue;
      }
      // angle cut
      //const double phi          = (mom_ion_phi < M_PI) ? mom_ion_phi : mom_ion_phi - 2 * M_PI;
      //const bool in_small_angle = (mom_ion_theta > m_thetaMinFullOMD && mom_ion_theta < m_thetaMaxFullOMD);
      //const bool in_large_angle = (mom_ion_theta > m_thetaMinPartialOMD && mom_ion_theta < m_thetaMaxPartialOMD);
      //if (!in_small_angle || (std::abs(phi) > 1 && !in_large_angle)) {
      //  continue;
      //}
      auto rc_part = smearMomentum(part);
      rc_part.status(kTagOMD);
      rc.push_back(rc_part);
      if (msgLevel(MSG::DEBUG)) {
        auto& rec_part = rc.back();
        const auto& part_p = part.getMomentum();
        const auto part_p_pt = std::hypot(part_p.x, part_p.y);
        const auto part_p_mag = std::hypot(part_p.x, part_p.y, part_p.z);
        const auto part_p_theta = std::atan2(std::hypot(part_p.x, part_p.y), part_p.z);
        debug() << fmt::format("Found OMD particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.getPDG(), part_p_mag, rec_part.momentum(), part_p_pt,
                               std::hypot(rec_part.p().x, rec_part.p().y), part_p_theta, rec_part.p().theta())
                << endmsg;
      }
    }
    return rc;
  }

  // all momentum smearing in EIC-smear for the far-forward region uses
  // the same 2 relations for P and Pt smearing (B0, RP, OMD)
  RecData smearMomentum(const edm4hep::ConstMCParticle& part) {
    const auto mom_ion = rotateLabToIonDirection(part.getMomentum());
    const double p     = std::hypot(mom_ion.x, mom_ion.y, mom_ion.z);
    const double dp    = (0.025 * p) * m_gaussDist();
    const double ps    = p + dp;

    const double pt  = std::hypot(mom_ion.x, mom_ion.y);
    const double dpt = (0.03 * pt) * m_gaussDist();
    // just apply relative smearing on px and py
    const double dpxs = (0.03 * mom_ion.x) * m_gaussDist(); //+ (1 + dpt / pt);
    const double dpys = (0.03 * mom_ion.y) * m_gaussDist(); //+ (1 + dpt / pt);

    const double pxs = mom_ion.x + dpxs; 
    const double pys = mom_ion.y + dpys; 

    // now get pz
    const double pzs = sqrt(ps * ps - pxs * pxs - pys * pys);

    // And build our 3-vector
    const edm4hep::Vector3f psmear_ion = {pxs, pys, pzs};
    const auto psmear = rotateIonToLabDirection(psmear_ion);
    const auto psmear_phi = std::atan2(psmear.x, psmear.y);
    const auto psmear_theta = std::atan2(std::hypot(psmear.x, psmear.y), psmear.z);
    eic::ReconstructedParticle rec_part;
    rec_part.ID({part.id(), algorithmID()});
    rec_part.p({psmear.x, psmear.y, psmear.z});
    rec_part.v({part.getVertex().x, part.getVertex().y, part.getVertex().z});
    rec_part.time(static_cast<float>(part.getTime()));
    rec_part.pid(part.getPDG());
    rec_part.status(0);
    rec_part.charge(static_cast<int16_t>(part.getCharge()));
    rec_part.weight(1.);
    rec_part.direction({psmear_theta, psmear_phi});
    rec_part.momentum(static_cast<float>(ps));
    rec_part.energy(std::hypot(ps, part.getMass()));
    rec_part.mass(static_cast<float>(part.getMass()));
    return rec_part;
  }

  // Rotate 25mrad about the y-axis
  edm4hep::Vector3f rotateLabToIonDirection(const edm4hep::Vector3f& vec) const {
    const auto sth = sin(-m_crossingAngle);
    const auto cth = cos(-m_crossingAngle);
    return {cth * vec.x + sth * vec.z, vec.y, -sth * vec.x + cth * vec.z};
  }

  edm4hep::Vector3f rotateIonToLabDirection(const edm4hep::Vector3f& vec) const {
    const auto sth = sin(m_crossingAngle);
    const auto cth = cos(m_crossingAngle);
    return {cth * vec.x + sth * vec.z, vec.y, -sth * vec.x + cth * vec.z};
  }

  edm4hep::Vector3f removeCrossingAngle(const edm4hep::Vector3f& vec) const {
    const auto sth = std::sin(-m_crossingAngle);
    const auto cth = std::cos(-m_crossingAngle);
    return {cth * vec.x + sth * vec.z, vec.y, -sth * vec.x + cth * vec.z};
  }

}; 

DECLARE_COMPONENT(SmearedFarForwardParticles)

} // namespace Jug::Fast

