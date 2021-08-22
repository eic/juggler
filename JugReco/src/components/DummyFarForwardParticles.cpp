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

namespace Jug::Reco {

class DummyFarForwardParticles : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  DataHandle<dd4pod::Geant4ParticleCollection> m_inputParticles{"inputMCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer,
                                                                     this};
  Gaudi::Property<bool> m_enableZDC{this, "enableZDC", true};
  Gaudi::Property<bool> m_enableB0{this, "enableB0", true};
  Gaudi::Property<bool> m_enableRP{this, "enableRP", true};
  Gaudi::Property<bool> m_enableOMD{this, "enableOMD", true};
  // RP default to 10-on-100 setting
  // P > 60% of beam energy (60% x 100GeV = 60GeV)
  // theta from 0.2mrad -> 5mrad
  Gaudi::Property<double> m_thetaMinRP{this, "thetaMinRP", 0.2e-3};
  Gaudi::Property<double> m_thetaMaxRP{this, "thetaMaxRP", 5e-3};
  Gaudi::Property<double> m_pMinRP{this, "pMinRP", 60};
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
  Gaudi::Property<double> m_pMinOMD{this, "pMinOMD", 25.};
  Gaudi::Property<double> m_pMaxOMD{this, "pMaxOMD", 60.};

  Rndm::Numbers m_gaussDist;

  DummyFarForwardParticles(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputCollection", m_inputParticles, "mcparticles");
    declareProperty("outputCollection", m_outputParticles, "ReconstructedParticles");
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

    std::vector<std::vector<eic::ReconstructedParticle>> rc_parts;
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
  std::vector<eic::ReconstructedParticle> zdc(const dd4pod::Geant4ParticleCollection& mc) {
    std::vector<eic::ReconstructedParticle> rc;
    for (const auto& part : mc) {
      if (part.genStatus() != 1) {
        continue;
      }
      // only detect neutrons and photons
      if (part.pdgID() != 2112 && part.pdgID() != 22) {
        continue;
      }
      // only 0-->4.5 mrad
      if (part.ps().theta() > 4.5 / 1000.) {
        continue;
      }
      const double E    = std::hypot(part.ps().mag(), part.mass());
      const double dE   = sqrt((0.05 * E) * (0.05 * E) + 0.5 * 0.5 * E) * m_gaussDist();
      const double Es   = E + dE;
      const double th   = part.ps().theta();
      const double dth  = (3e-3 / sqrt(E)) * m_gaussDist();
      const double ths  = th + dth;
      const double phi  = part.ps().phi();
      const double dphi = 0;
      const double phis = phi + dphi;

      const double moms = sqrt(Es * Es - part.mass() * part.mass());
      const eic::VectorPolar mom3s{moms, ths, phis};
      eic::ReconstructedParticle rec_part{part.ID(),
                                          mom3s,
                                          {part.vs().x, part.vs().y, part.vs().z},
                                          static_cast<float>(part.time()),
                                          part.pdgID(),
                                          0,
                                          static_cast<int16_t>(part.charge()),
                                          algorithmID(),
                                          1.,
                                          static_cast<float>(moms),
                                          static_cast<float>(Es),
                                          static_cast<float>(part.mass())};
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
  std::vector<eic::ReconstructedParticle> b0(const dd4pod::Geant4ParticleCollection& mc) {
    std::vector<eic::ReconstructedParticle> rc;
    for (const auto& part : mc) {
      if (part.genStatus() != 1) {
        continue;
      }
      // only detect charged hadrons and photons
      if (part.pdgID() != 2212 && part.pdgID() != -2212 && part.pdgID() != 211 && part.pdgID() != -211 &&
          part.pdgID() != 321 && part.pdgID() != -321 && part.pdgID() != 22) {
        continue;
      }
      // only 6-->20 mrad
      if (part.ps().theta() < m_thetaMinB0 || part.ps().theta() > m_thetaMaxB0) {
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

  std::vector<eic::ReconstructedParticle> rp(const dd4pod::Geant4ParticleCollection& mc) {
    std::vector<eic::ReconstructedParticle> rc;
    for (const auto& part : mc) {
      if (part.genStatus() != 1) {
        continue;
      }
      // only detect protons
      if (part.pdgID() != 2212) {
        continue;
      }
      if (part.ps().theta() < m_thetaMinRP || part.ps().theta() > m_thetaMaxRP || part.ps().mag() < m_pMinRP) {
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

  std::vector<eic::ReconstructedParticle> omd(const dd4pod::Geant4ParticleCollection& mc) {
    std::vector<eic::ReconstructedParticle> rc;
    for (const auto& part : mc) {
      if (part.genStatus() != 1) {
        continue;
      }
      // only detect protons
      if (part.pdgID() != 2212) {
        continue;
      }
      // momentum cut
      if (part.ps().mag() < m_pMinOMD || part.ps().mag() > m_pMaxOMD) {
        continue;
      }
      // angle cut
      const double phi          = (part.ps().phi() < M_PI) ? part.ps().phi() : part.ps().phi() - 2 * M_PI;
      const bool in_small_angle = (part.ps().theta() > m_thetaMinFullOMD && part.ps().theta() < m_thetaMaxFullOMD);
      const bool in_large_angle =
          (part.ps().theta() > m_thetaMinPartialOMD && part.ps().theta() < m_thetaMaxPartialOMD);
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
  eic::ReconstructedParticle smearMomentum(const dd4pod::ConstGeant4Particle& part) {
    const double p  = part.ps().mag();
    const double dp = (0.005 * p) * m_gaussDist();
    const double ps = p + dp;

    const double pt  = std::hypot(part.ps().x, part.ps().y);
    const double dpt = (0.03 * pt) * m_gaussDist();
    // just apply relative smearing on px and py
    const double pxs = part.ps().x + (1 + dpt / pt);
    const double pys = part.ps().y + (1 + dpt / pt);
    // now get pz
    const double pzs      = sqrt(ps * ps - pxs * pxs - pys * pys);
    eic::VectorXYZ psmear = {pxs, pys, pzs};
    return {part.ID(),
            psmear,
            {part.vs().x, part.vs().y, part.vs().z},
            static_cast<float>(part.time()),
            part.pdgID(),
            0,
            static_cast<int16_t>(part.charge()),
            algorithmID(),
            1.,
            static_cast<float>(ps),
            static_cast<float>(std::hypot(psmear.mag(), part.mass())),
            static_cast<float>(part.mass())};
  }

}; // namespace Jug::Reco

DECLARE_COMPONENT(DummyFarForwardParticles)

} // namespace Jug::Reco

