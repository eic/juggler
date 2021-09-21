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
#include "eicd/ReconstructedParticleRelationsCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/VectorPolar.h"

namespace Jug::Reco {

class ParticlesWithTruthPID : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  DataHandle<dd4pod::Geant4ParticleCollection> m_inputTruthCollection{"inputMCParticles", Gaudi::DataHandle::Reader,
                                                                      this};
  DataHandle<eic::TrackParametersCollection> m_inputTrackCollection{"inputTrackParameters", Gaudi::DataHandle::Reader,
                                                                    this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticleCollection{
      "ReconstructedParticles", Gaudi::DataHandle::Writer, this};
  DataHandle<eic::ReconstructedParticleRelationsCollection> m_outputRelationsCollection{
      "ReconstructedParticleRelations", Gaudi::DataHandle::Writer, this};

  // Matching momentum tolerance requires 10% by default;
  Gaudi::Property<double> m_pRelativeTolerance{this, "pRelativeTolerance", {0.1}};
  // Matching phi tolerance of 10 mrad
  Gaudi::Property<double> m_phiTolerance{this, "phiTolerance", {0.030}};
  // Matchin eta tolerance of 0.1
  Gaudi::Property<double> m_etaTolerance{this, "etaTolerance", {0.2}};

  const int32_t m_kMonteCarloSource{uniqueID<int32_t>("mcparticles")};

  ParticlesWithTruthPID(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputTruthCollection, "mcparticles");
    declareProperty("inputTrackParameters", m_inputTrackCollection, "outputTrackParameters");
    declareProperty("outputParticles", m_outputParticleCollection, "ReconstructedParticles");
    declareProperty("outputRelations", m_outputRelationsCollection, "ReconstructedParticleRelations");
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
    auto& rel          = *(m_outputRelationsCollection.createAndPut());

    const double sinPhiOver2Tolerance = sin(0.5 * m_phiTolerance);
    std::vector<bool> consumed(mc.size(), false);
    int ID = 0;
    for (const auto& trk : tracks) {
      const eic::VectorXYZ mom =
          eic::VectorPolar(1.0 / std::abs(trk.qOverP()), trk.direction().theta, trk.direction().phi);
      const auto charge_rec = std::copysign(1., trk.qOverP());
      // utility variables for matching
      int best_match    = -1;
      double best_delta = std::numeric_limits<double>::max();
      for (size_t ip = 0; ip < mc.size(); ++ip) {
        const auto& mcpart = mc[ip];
        if (consumed[ip] || mcpart.genStatus() > 1 || mcpart.charge() == 0 || mcpart.charge() * charge_rec < 0) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "ignoring non-primary/neutral/opposite charge particle" << endmsg;
          }
          continue;
        }
        const double dp_rel = std::abs((mom.mag() - mcpart.ps().mag()) / mcpart.ps().mag());
        // check the tolerance for sin(dphi/2) to avoid the hemisphere problem and allow
        // for phi rollovers
        const double dsphi = std::abs(sin(0.5 * (mom.phi() - mcpart.ps().phi())));
        const double deta  = std::abs((mom.eta() - mcpart.ps().eta()));

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
        best_pid             = mcpart.pdgID();
        vertex               = {mcpart.vs().x, mcpart.vs().y, mcpart.vs().z};
        time                 = mcpart.time();
        mass                 = mcpart.mass();
        mcID                 = {mcpart.ID(), m_kMonteCarloSource};
      }
      eic::ReconstructedParticle rec_part{{ID++, algorithmID()},
                                          mom,
                                          vertex,
                                          time,
                                          best_pid,
                                          static_cast<int16_t>(best_match >= 0 ? 0 : -1) /* status */,
                                          static_cast<int16_t>(charge_rec),
                                          1. /* weight */,
                                          {mom.theta(), mom.phi()},
                                          mom.mag(),
                                          std::hypot(mom.mag(), mass),
                                          mass};
      eic::ReconstructedParticleRelations rel_part;
      rel_part.recID(rec_part.ID());
      rel_part.mcID(mcID);
      part.push_back(rec_part);
      rel.push_back(rel_part);
      if (msgLevel(MSG::DEBUG)) {
        if (best_match > 0) {
          const auto& mcpart = mc[best_match];
          debug() << fmt::format("Matched track {} with MC particle {}\n", trk.ID().value, best_match) << endmsg;
          debug() << fmt::format("  - Track: (mom: {}, theta: {}, phi: {}, charge: {})", mom.mag(), mom.theta(),
                                 mom.phi(), charge_rec)
                  << endmsg;
          debug() << fmt::format("  - MC particle: (mom: {}, theta: {}, phi: {}, charge: {}, type: {}",
                                 mcpart.ps().mag(), mcpart.ps().theta(), mcpart.ps().phi(), mcpart.charge(),
                                 mcpart.pdgID())
                  << endmsg;
        } else {
          debug() << fmt::format("Did not find a good match for track {} \n", trk.ID().value) << endmsg;
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

} // namespace Jug::Reco

