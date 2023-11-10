// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

#include <algorithm>
#include <cmath>
#include <fmt/format.h>

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4hep/utils/vector_utils.h"

namespace {
enum DetectorTags { kTagB0 = 1, kTagRP = 2, kTagOMD = 3, kTagZDC = 4 };
}

namespace Jug::Fast {

class SmearedFarForwardParticles : public GaudiAlgorithm {
private:
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"inputMCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::ReconstructedParticleCollection> m_outputParticles{"SmearedFarForwardParticles",
                                                                      Gaudi::DataHandle::Writer, this};
  DataHandle<edm4eic::MCRecoParticleAssociationCollection> m_outputAssocCollection{"MCRecoParticleAssociation",
                                                                                Gaudi::DataHandle::Writer, this};

  Gaudi::Property<bool> m_enableZDC{this, "enableZDC", true};
  Gaudi::Property<bool> m_enableB0{this, "enableB0", true};
  Gaudi::Property<bool> m_enableRP{this, "enableRP", true};
  Gaudi::Property<bool> m_enableOMD{this, "enableOMD", true};

  // Beam energy, only used to determine the RP/OMD momentum ranges
  Gaudi::Property<double> m_ionBeamEnergy{this, "ionBeamEnergy", 0.};
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
  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle",
                                          -0.025}; //-0.025}; -- causes double rotation with afterburner

  Rndm::Numbers m_gaussDist;

  using RecPart = edm4eic::MutableReconstructedParticle;
  using Assoc   = edm4eic::MutableMCRecoParticleAssociation;
  using RecData = std::pair<RecPart, Assoc>;

public:
  SmearedFarForwardParticles(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticles, "MCParticles");
    declareProperty("outputParticles", m_outputParticles, "ReconstructedParticles");
    declareProperty("outputAssociations", m_outputAssocCollection, "MCRecoParticleAssociation");
  }
  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
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
    auto& assoc    = *(m_outputAssocCollection.createAndPut());

    double ionBeamEnergy = 0;
    if (m_ionBeamEnergy > 0) {
      ionBeamEnergy = m_ionBeamEnergy;
    } else {
      for (const auto& part : mc) {
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
      for (const auto& [part, link] : det) {
        rc.push_back(part);
        assoc.push_back(link);
      }
    }
    return StatusCode::SUCCESS;
  }

private:
  // ZDC smearing as in eic_smear
  // https://github.com/eic/eicsmeardetectors/blob/9a1831dd97bf517b80a06043b9ee4bfb96b483d8/SmearMatrixDetector_0_1_FF.cxx#L224
  std::vector<RecData> zdc(const edm4hep::MCParticleCollection& mc, const double /* ionBeamEnergy */) {
    std::vector<RecData> rc;
    for (const auto& part : mc) {
      if (part.getGeneratorStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with generatorStatus = " << part.getGeneratorStatus() << endmsg;
        }
        continue;
      }
      // only detect neutrons and photons
      const auto mom_ion = rotateLabToIonDirection(part.getMomentum());
      if (part.getPDG() != 2112 && part.getPDG() != 22) {
        continue;
      }
      // only 0-->4.5 mrad
      const double mom_ion_theta = edm4hep::utils::anglePolar(mom_ion);
      const double mom_ion_phi   = edm4hep::utils::angleAzimuthal(mom_ion);
      if (mom_ion_theta > 4.5 / 1000.) {
        continue;
      }

      double conTerm = 0.05;  // default 5%
      double stoTerm = 0.5;   // default 50%
      double angTerm = 0.003; // 3mrad

      if (part.getPDG() == 2112) {
        conTerm = 0.05;                 // default 5%
        stoTerm = 0.5;                  // default 50%
        angTerm = 0.003;                // 3mrad
      } else if (part.getPDG() == 22) { // EMCAL expected to have slightly better performance
        conTerm = 0.03;                 // default 3%
        stoTerm = 0.10;                 // default 10% for WSciFi
        angTerm = 0.001;                // 1mrad is the detault for the block size
      }

      // explicit double precision due to E*E - m*m
      const double E    = part.getEnergy();
      const double dE   = sqrt((conTerm * E) * (conTerm * E) + stoTerm * stoTerm * E) * m_gaussDist(); // 50%/SqrtE + 5%
      const double Es   = E + dE;
      const double th   = mom_ion_theta;
      const double dth  = (angTerm / sqrt(E)) * m_gaussDist();
      const double ths  = th + dth;
      const double phi  = mom_ion_phi;
      const double dphi = 0;
      const double phis = phi + dphi;
      const double moms = sqrt(Es * Es - part.getMass() * part.getMass());
      // now cast back into float
      const auto mom3s_ion = edm4hep::utils::sphericalToVector(moms, ths, phis);
      const auto mom3s     = rotateIonToLabDirection(mom3s_ion);
      RecPart rec_part;
      rec_part.setType(kTagZDC);
      rec_part.setEnergy(static_cast<float>(Es));
      rec_part.setMomentum({mom3s.x, mom3s.y, mom3s.z});
      rec_part.setReferencePoint({static_cast<float>(part.getVertex().x), static_cast<float>(part.getVertex().y),
                                  static_cast<float>(part.getVertex().z)});
      rec_part.setCharge(static_cast<int16_t>(part.getCharge()));
      rec_part.setMass(static_cast<float>(part.getMass()));
      rec_part.setGoodnessOfPID(1.);
      rec_part.setPDG(part.getPDG());
      Assoc assoc;
      assoc.setRecID(rec_part.getObjectID().index);
      assoc.setSimID(part.getObjectID().index);
      assoc.setWeight(1.);
      assoc.setRec(rec_part);
      //assoc.setSim(part);

      // rec_part.mcID();
      rc.emplace_back(rec_part, assoc);

      if (msgLevel(MSG::DEBUG)) {
        const auto& part_p    = part.getMomentum();
        const auto part_p_mag = std::hypot(part_p.x, part_p.y, part_p.z);
        debug()
            << fmt::format(
                   "Found ZDC particle: {}, Etrue: {}, Emeas: {}, ptrue: {}, pmeas: {}, theta_true: {}, theta_meas: {}",
                   part.getPDG(), E, rec_part.getEnergy(), part_p_mag, edm4hep::utils::magnitude(rec_part.getMomentum()), th,
                   edm4hep::utils::anglePolar(rec_part.getMomentum()))
            << endmsg;
      }
    }
    return rc;
  }
  // Fast B0 as in
  // https://github.com/eic/eicsmeardetectors/blob/9a1831dd97bf517b80a06043b9ee4bfb96b483d8/SmearMatrixDetector_0_1_FF.cxx#L254
  std::vector<RecData> b0(const edm4hep::MCParticleCollection& mc, const double /* ionBeamEnergy */) {
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
      const auto mom_ion = removeCrossingAngle(part.getMomentum()); // rotateLabToIonDirection(part.getMomentum());
      const auto mom_ion_theta = edm4hep::utils::anglePolar(mom_ion);
      if (mom_ion_theta < m_thetaMinB0 || mom_ion_theta > m_thetaMaxB0) {
        continue;
      }
      auto [rc_part, assoc] = smearMomentum(part);
      // we don't detect photon energy, just its angles and presence
      if (part.getPDG() == 22) {
        rc_part.setMomentum({0, 0, 0});
        rc_part.setEnergy(0);
      }
      rc_part.setType(kTagB0);
      rc.emplace_back(rc_part, assoc);
      if (msgLevel(MSG::DEBUG)) {
        const auto& part_p      = part.getMomentum();
        const auto part_p_pt    = edm4hep::utils::magnitudeTransverse(part_p);
        const auto part_p_mag   = edm4hep::utils::magnitude(part_p);
        const auto part_p_theta = edm4hep::utils::anglePolar(part_p);
        debug() << fmt::format("Found B0 particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.getPDG(), part_p_mag, edm4hep::utils::magnitude(rc_part.momentum()), part_p_pt,
                               edm4hep::utils::magnitudeTransverse(rc_part.momentum()), part_p_theta,
                               edm4hep::utils::anglePolar(rc_part.momentum()))
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
      const auto mom_ion = removeCrossingAngle(part.getMomentum()); // rotateLabToIonDirection(part.getMomentum());
      const auto mom_ion_theta = edm4hep::utils::anglePolar(mom_ion);
      if (mom_ion_theta < m_thetaMinRP || mom_ion_theta > m_thetaMaxRP ||
          mom_ion.z < m_pMinRigidityRP * ionBeamEnergy) {
        continue;
      }
      auto [rc_part, assoc] = smearMomentum(part);
      rc_part.setType(kTagRP);
      rc.emplace_back(rc_part, assoc);
      if (msgLevel(MSG::DEBUG)) {
        const auto& part_p      = part.getMomentum();
        const auto part_p_pt    = edm4hep::utils::magnitudeTransverse(part_p);
        const auto part_p_mag   = edm4hep::utils::magnitude(part_p);
        const auto part_p_theta = edm4hep::utils::anglePolar(part_p);
        debug() << fmt::format("Found RP particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.getPDG(), part_p_mag, edm4hep::utils::magnitude(rc_part.momentum()), part_p_pt,
                               edm4hep::utils::magnitudeTransverse(rc_part.momentum()), part_p_theta,
                               edm4hep::utils::anglePolar(rc_part.momentum()))
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
      const auto mom_ion = removeCrossingAngle(part.getMomentum()); // rotateLabToIonDirection(part.getMomentum());
      if (mom_ion.z < m_pMinRigidityOMD * ionBeamEnergy || mom_ion.z > m_pMaxRigidityOMD * ionBeamEnergy) {
        continue;
      }
      auto [rc_part, assoc] = smearMomentum(part);
      rc_part.setType(kTagOMD);
      rc.emplace_back(rc_part, assoc);
      if (msgLevel(MSG::DEBUG)) {
        const auto& part_p      = part.getMomentum();
        const auto part_p_pt    = edm4hep::utils::magnitudeTransverse(part_p);
        const auto part_p_mag   = edm4hep::utils::magnitude(part_p);
        const auto part_p_theta = edm4hep::utils::anglePolar(part_p);
        debug() << fmt::format("Found OMD particle: {}, ptrue: {}, pmeas: {}, pttrue: {}, ptmeas: {}, theta_true: {}, "
                               "theta_meas: {}",
                               part.getPDG(), part_p_mag, edm4hep::utils::magnitude(rc_part.momentum()), part_p_pt,
                               edm4hep::utils::magnitudeTransverse(rc_part.momentum()), part_p_theta,
                               edm4hep::utils::anglePolar(rc_part.momentum()))
                << endmsg;
      }
    }
    return rc;
  }

  // all momentum smearing in EIC-smear for the far-forward region uses
  // the same 2 relations for P and Pt smearing (B0, RP, OMD)
  RecData smearMomentum(const edm4hep::MCParticle& part) {
    const auto mom_ion = rotateLabToIonDirection(part.getMomentum());
    const double p     = std::hypot(mom_ion.x, mom_ion.y, mom_ion.z);
    const double dp    = (0.025 * p) * m_gaussDist();
    const double ps    = p + dp;

    // const double pt  = std::hypot(mom_ion.x, mom_ion.y);
    // const double dpt = (0.03 * pt) * m_gaussDist();
    // just apply relative smearing on px and py
    const double dpxs = (0.03 * mom_ion.x) * m_gaussDist(); //+ (1 + dpt / pt);
    const double dpys = (0.03 * mom_ion.y) * m_gaussDist(); //+ (1 + dpt / pt);

    const double pxs = mom_ion.x + dpxs;
    const double pys = mom_ion.y + dpys;

    // now get pz
    const double pzs = sqrt(ps * ps - pxs * pxs - pys * pys);

    // And build our 3-vector
    const edm4hep::Vector3f psmear_ion{static_cast<float>(pxs), static_cast<float>(pys), static_cast<float>(pzs)};
    const auto psmear = rotateIonToLabDirection(psmear_ion);
    edm4eic::MutableReconstructedParticle rec_part;
    rec_part.setType(-1);
    rec_part.setEnergy(std::hypot(ps, part.getMass()));
    rec_part.setMomentum({psmear.x, psmear.y, psmear.z});
    rec_part.setReferencePoint({static_cast<float>(part.getVertex().x), static_cast<float>(part.getVertex().y),
                                static_cast<float>(part.getVertex().z)});
    rec_part.setCharge(static_cast<int16_t>(part.getCharge()));
    rec_part.setMass(static_cast<float>(part.getMass()));
    rec_part.setGoodnessOfPID(1); // perfect PID
    rec_part.setPDG(part.getPDG());
    Assoc assoc;
    assoc.setRecID(rec_part.getObjectID().index);
    assoc.setSimID(part.getObjectID().index);
    assoc.setWeight(1.);
    assoc.setRec(rec_part);
    //assoc.setSim(part);

    return {rec_part, assoc};
  }

  // Rotate 25mrad about the y-axis
  edm4hep::Vector3f rotateLabToIonDirection(const edm4hep::Vector3f& vec) const {
    const auto sth = sin(-m_crossingAngle);
    const auto cth = cos(-m_crossingAngle);
    return {static_cast<float>(cth * vec.x + sth * vec.z), static_cast<float>(vec.y),
            static_cast<float>(-sth * vec.x + cth * vec.z)};
  }

  edm4hep::Vector3f rotateIonToLabDirection(const edm4hep::Vector3f& vec) const {
    const auto sth = sin(m_crossingAngle);
    const auto cth = cos(m_crossingAngle);
    return {static_cast<float>(cth * vec.x + sth * vec.z), static_cast<float>(vec.y),
            static_cast<float>(-sth * vec.x + cth * vec.z)};
  }

  edm4hep::Vector3f removeCrossingAngle(const edm4hep::Vector3f& vec) const {
    const auto sth = std::sin(-m_crossingAngle);
    const auto cth = std::cos(-m_crossingAngle);
    return {static_cast<float>(cth * vec.x + sth * vec.z), static_cast<float>(vec.y),
            static_cast<float>(-sth * vec.x + cth * vec.z)};
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(SmearedFarForwardParticles)

} // namespace Jug::Fast
