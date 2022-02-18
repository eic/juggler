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
#include "eicd/TrackParametersCollection.h"
#include "eicd/VectorPolar.h"

namespace Jug::Fast {

class ParticlesWithTruthPID : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  DataHandle<edm4hep::MCParticleCollection> m_inputTruthCollection{"inputMCParticles", Gaudi::DataHandle::Reader,
                                                                      this};
  DataHandle<eic::TrackParametersCollection> m_inputTrackCollection{"inputTrackParameters", Gaudi::DataHandle::Reader,
                                                                    this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticleCollection{
      "ReconstructedParticles", Gaudi::DataHandle::Writer, this};

  // Matching momentum tolerance requires 10% by default;
  Gaudi::Property<double> m_pRelativeTolerance{this, "pRelativeTolerance", {0.1}};
  // Matching phi tolerance of 10 mrad
  Gaudi::Property<double> m_phiTolerance{this, "phiTolerance", {0.030}};
  // Matchin eta tolerance of 0.1
  Gaudi::Property<double> m_etaTolerance{this, "etaTolerance", {0.2}};

  const int32_t m_kMonteCarloSource{uniqueID<int32_t>("MCParticles")};

  ParticlesWithTruthPID(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputTruthCollection, "MCParticles");
    declareProperty("inputTrackParameters", m_inputTrackCollection, "outputTrackParameters");
    declareProperty("outputParticles", m_outputParticleCollection, "ReconstructedParticles");
  }
  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;

    return StatusCode::SUCCESS;
  }
  StatusCode execute() override {
    // input collection
    const auto& mc     = *(m_inputTruthCollection.get());
    const auto& tracks = *(m_inputTrackCollection.get());
    auto& part         = *(m_outputParticleCollection.createAndPut());

    const double sinPhiOver2Tolerance = sin(0.5 * m_phiTolerance);
    std::vector<bool> consumed(mc.size(), false);
    int ID = 0;
    for (const auto& trk : tracks) {
      const eic::VectorXYZ mom =
          eic::VectorPolar(1.0 / std::abs(trk.qOverP()), trk.theta(), trk.phi());
      const auto charge_rec = std::copysign(1., trk.qOverP());
      // utility variables for matching
      int best_match    = -1;
      double best_delta = std::numeric_limits<double>::max();
      for (size_t ip = 0; ip < mc.size(); ++ip) {
        const auto& mcpart = mc[ip];
        if (consumed[ip] || mcpart.getGeneratorStatus() > 1 || mcpart.getCharge() == 0 || mcpart.getCharge() * charge_rec < 0) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "ignoring non-primary/neutral/opposite charge particle" << endmsg;
          }
          continue;
        }
        const auto& p = mcpart.getMomentum();
        const auto p_mag = std::hypot(p.x, p.y, p.z);
        const auto p_phi = std::atan2(p.y, p.x);
        const auto p_eta = std::atanh(p.z / p_mag);
        const double dp_rel = std::abs((mom.mag() - p_mag) / p_mag);
        // check the tolerance for sin(dphi/2) to avoid the hemisphere problem and allow
        // for phi rollovers
        const double dsphi = std::abs(sin(0.5 * (mom.phi() - p_phi)));
        const double deta  = std::abs((mom.eta() - p_eta));

        if (dp_rel < m_pRelativeTolerance && deta < m_etaTolerance && dsphi < sinPhiOver2Tolerance) {
          const double delta =
              std::hypot(dp_rel / m_pRelativeTolerance, deta / m_etaTolerance, dsphi / sinPhiOver2Tolerance);
          if (delta < best_delta) {
            best_match = ip;
            best_delta = delta;
          }
        }
      }
      int32_t best_pid = 0;
      eic::VectorXYZ vertex;
      float time = 0;
      float mass = 0;
      eic::Index mcID;
      if (best_match >= 0) {
        consumed[best_match] = true;
        const auto& mcpart   = mc[best_match];
        best_pid             = mcpart.getPDG();
        vertex               = {mcpart.getVertex().x, mcpart.getVertex().y, mcpart.getVertex().z};
        time                 = mcpart.getTime();
        mass                 = mcpart.getMass();
      }
      auto rec_part = part.create();
      rec_part.ID({ID++, algorithmID()});
      rec_part.p(mom);
      rec_part.v(vertex);
      rec_part.time(time);
      rec_part.pid(best_pid);
      rec_part.status(static_cast<int16_t>(best_match >= 0 ? 0 : -1));
      rec_part.charge(static_cast<int16_t>(charge_rec));
      rec_part.weight(1.);
      rec_part.direction({mom.theta(), mom.phi()});
      rec_part.momentum(mom.mag());
      rec_part.energy(std::hypot(mom.mag(), mass));
      rec_part.mass(mass);
      rec_part.mcID(mcID);

      if (msgLevel(MSG::DEBUG)) {
        if (best_match > 0) {
          const auto& mcpart = mc[best_match];
          debug() << fmt::format("Matched track {} with MC particle {}\n", ID, best_match) << endmsg;
          debug() << fmt::format("  - Track: (mom: {}, theta: {}, phi: {}, charge: {})", mom.mag(), mom.theta(),
                                 mom.phi(), charge_rec)
                  << endmsg;
          const auto& p = mcpart.getMomentum();
          const auto p_mag = std::hypot(p.x, p.y, p.z);
          const auto p_phi = std::atan2(p.y, p.x);
          const auto p_theta = std::atan2(std::hypot(p.x, p.y), p.z);
          debug() << fmt::format("  - MC particle: (mom: {}, theta: {}, phi: {}, charge: {}, type: {}",
                                 p_mag, p_theta, p_phi, mcpart.getCharge(), mcpart.getPDG())
                  << endmsg;
        } else {
          debug() << fmt::format("Did not find a good match for track {} \n", ID) << endmsg;
          debug() << fmt::format("  - Track: (mom: {}, theta: {}, phi: {}, charge: {})", mom.mag(), mom.theta(),
                                 mom.phi(), charge_rec)
                  << endmsg;
        }
      }
    }

    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(ParticlesWithTruthPID)

} // namespace Jug::Fast

