#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/VectorPolar.h"

namespace Jug::Rec {

class ParticleWithTruthPID : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  DataHandle<dd4pod::Geant4ParticleCollection> m_inputTruthCollection{"inputMC", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::TrackParametersCollection> m_inputTrackCollection{"inputTracks", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticleCollection{"outputParticles",
                                                                              Gaudi::DataHandle::Writer, this};
  ParticleWithTruthPID(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMC", m_inputTruthCollection, "mcparticles");
    declareProperty("inputTracks", m_inputTrackCollection, "outputTrackParameters");
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

    std::vector<bool> consumed(tracks.size(), false);
    int ID = 0;
    for (const auto& p : mc) {
      if (p.genStatus() != 1) {
        continue;
      }
      // check if we find a good match
      int best_match    = -1;
      double best_delta = std::numeric_limits<double>::max();
      for (int it = 0; it < tracks.size(); ++it) {
        if (consumed[it]) {
          continue;
        }
        const auto& trk = tracks[it];
        eic::VectorPolar mom{1.0 / std::abs(trk.qOverP()), trk.direction().theta, trk.direction().phi};
        double delta = std::hypot(p.ps().x - mom.x(), p.ps().y - mom.y(), p.ps().z - mom.z());
        if (delta < best_delta) {
          best_match = it;
          best_delta = delta;
        }
      }
      if (best_match > 0) {
        consumed[best_match] = true;
        const auto& trk      = tracks[best_match];
        eic::VectorPolar mom{1.0 / std::abs(trk.qOverP()), trk.direction().theta, trk.direction().phi};
        eic::ReconstructedParticle rec_part{ID++,                                                  // index
                                            {mom.x(), mom.y(), mom.z()},                           // momentum
                                            {0., 0., 0.},                                          // vertex
                                            0.,                                                    // time
                                            p.pdgID(),                                             // PID
                                            static_cast<int16_t>(0),                               // Status
                                            static_cast<int16_t>(std::copysign(1., trk.qOverP())), // charge
                                            algorithmID(),                                         // Algorithm type
                                            1.,                                                    // particle weight
                                            mom.mag(),                       // 3-momentum magnitude [GeV]
                                            std::hypot(mom.mag(), p.mass()), // energy [GeV]
                                            static_cast<float>(p.mass())};   // mass [GeV]
        part.push_back(rec_part);
      }
    }

    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(ParticleWithTruthPID)

} // namespace Jug::Rec

