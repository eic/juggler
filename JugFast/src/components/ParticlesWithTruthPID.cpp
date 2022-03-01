#include <algorithm>
#include <cmath>

#include <fmt/format.h>

#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/MCRecoParticleAssociationCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/vector_utils.h"

namespace Jug::Fast {

class ParticlesWithTruthPID : public GaudiAlgorithm {
public:
  DataHandle<edm4hep::MCParticleCollection> m_inputTruthCollection{"inputMCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::TrackParametersCollection> m_inputTrackCollection{"inputTrackParameters", Gaudi::DataHandle::Reader,
                                                                     this};
  DataHandle<eicd::ReconstructedParticleCollection> m_outputParticleCollection{"ReconstructedParticles",
                                                                               Gaudi::DataHandle::Writer, this};
  DataHandle<eicd::MCRecoParticleAssociationCollection> m_outputAssocCollection{"MCRecoParticleAssociation",
                                                                                Gaudi::DataHandle::Writer, this};

  // Matching momentum tolerance requires 10% by default;
  Gaudi::Property<double> m_pRelativeTolerance{this, "pRelativeTolerance", {0.1}};
  // Matching phi tolerance of 10 mrad
  Gaudi::Property<double> m_phiTolerance{this, "phiTolerance", {0.030}};
  // Matchin eta tolerance of 0.1
  Gaudi::Property<double> m_etaTolerance{this, "etaTolerance", {0.2}};

  ParticlesWithTruthPID(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputTruthCollection, "MCParticles");
    declareProperty("inputTrackParameters", m_inputTrackCollection, "outputTrackParameters");
    declareProperty("outputParticles", m_outputParticleCollection, "ReconstructedParticles");
    declareProperty("outputAssociations", m_outputAssocCollection, "MCRecoParticleAssociation");
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
    auto& assoc        = *(m_outputAssocCollection.createAndPut());

    const double sinPhiOver2Tolerance = sin(0.5 * m_phiTolerance);
    std::vector<bool> consumed(mc.size(), false);
    for (const auto& trk : tracks) {
      const auto mom        = eicd::sphericalToVector(1.0 / std::abs(trk.getQOverP()), trk.getTheta(), trk.getPhi());
      const auto charge_rec = trk.getCharge();
      // utility variables for matching
      int best_match    = -1;
      double best_delta = std::numeric_limits<double>::max();
      for (size_t ip = 0; ip < mc.size(); ++ip) {
        const auto& mcpart = mc[ip];
        if (consumed[ip] || mcpart.getGeneratorStatus() > 1 || mcpart.getCharge() == 0 ||
            mcpart.getCharge() * charge_rec < 0) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "ignoring non-primary/neutral/opposite charge particle" << endmsg;
          }
          continue;
        }
        const auto& p       = mcpart.getMomentum();
        const auto p_mag    = std::hypot(p.x, p.y, p.z);
        const auto p_phi    = std::atan2(p.y, p.x);
        const auto p_eta    = std::atanh(p.z / p_mag);
        const double dp_rel = std::abs((eicd::magnitude(mom) - p_mag) / p_mag);
        // check the tolerance for sin(dphi/2) to avoid the hemisphere problem and allow
        // for phi rollovers
        const double dsphi = std::abs(sin(0.5 * (eicd::angleAzimuthal(mom) - p_phi)));
        const double deta  = std::abs((eicd::eta(mom) - p_eta));

        if (dp_rel < m_pRelativeTolerance && deta < m_etaTolerance && dsphi < sinPhiOver2Tolerance) {
          const double delta =
              std::hypot(dp_rel / m_pRelativeTolerance, deta / m_etaTolerance, dsphi / sinPhiOver2Tolerance);
          if (delta < best_delta) {
            best_match = ip;
            best_delta = delta;
          }
        }
      }
      auto rec_part       = part.create();
      int32_t best_pid    = 0;
      auto referencePoint = rec_part.referencePoint();
      // float time          = 0;
      float mass = 0;
      if (best_match >= 0) {
        consumed[best_match] = true;
        const auto& mcpart   = mc[best_match];
        best_pid             = mcpart.getPDG();
        referencePoint       = {
            static_cast<float>(mcpart.getVertex().x), static_cast<float>(mcpart.getVertex().y),
            static_cast<float>(mcpart.getVertex().z)}; // @TODO: not sure if vertex/reference poitn makes sense here
        // time                 = mcpart.getTime();
        mass = mcpart.getMass();
      }
      rec_part.setType(static_cast<int16_t>(best_match >= 0 ? 0 : -1)); // @TODO: determine type codes
      rec_part.setEnergy(std::hypot(eicd::magnitude(mom), mass));
      rec_part.setMomentum(mom);
      rec_part.setReferencePoint(referencePoint);
      rec_part.setCharge(charge_rec);
      rec_part.setMass(mass);
      rec_part.setGoodnessOfPID(1); // perfect PID
      rec_part.setPDG(best_pid);
      // rec_part.covMatrix()  // @TODO: covariance matrix on 4-momentum
      // Also write MC <--> truth particle association
      auto rec_assoc = assoc.create();
      rec_assoc.setSimID(best_match);
      rec_assoc.setWeight(1);
      rec_assoc.setRec(rec_part);

      if (msgLevel(MSG::DEBUG)) {
        if (best_match > 0) {
          const auto& mcpart = mc[best_match];
          debug() << fmt::format("Matched track with MC particle {}\n", best_match) << endmsg;
          debug() << fmt::format("  - Track: (mom: {}, theta: {}, phi: {}, charge: {})", eicd::magnitude(mom),
                                 eicd::anglePolar(mom), eicd::angleAzimuthal(mom), charge_rec)
                  << endmsg;
          const auto& p      = mcpart.getMomentum();
          const auto p_mag   = eicd::magnitude(p);
          const auto p_phi   = eicd::angleAzimuthal(p);
          const auto p_theta = eicd::anglePolar(p);
          debug() << fmt::format("  - MC particle: (mom: {}, theta: {}, phi: {}, charge: {}, type: {}", p_mag, p_theta,
                                 p_phi, mcpart.getCharge(), mcpart.getPDG())
                  << endmsg;
        } else {
          debug() << fmt::format("Did not find a good match for track \n") << endmsg;
          debug() << fmt::format("  - Track: (mom: {}, theta: {}, phi: {}, charge: {})", eicd::magnitude(mom),
                                 eicd::anglePolar(mom), eicd::angleAzimuthal(mom), charge_rec)
                  << endmsg;
        }
      }
    }

    return StatusCode::SUCCESS;
  }
}; // namespace Jug::Fast

DECLARE_COMPONENT(ParticlesWithTruthPID)

} // namespace Jug::Fast

