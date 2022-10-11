// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten

#include <cmath>
// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugTrack/Track.hpp"

#include "edm4eic/ClusterCollection.h"
#include "edm4eic/TrackerHitCollection.h"
#include "edm4eic/vector_utils.h"

#include "Acts/Surfaces/PerigeeSurface.hpp"

///// (Reconstructed) track parameters e.g. close to the vertex.
// using TrackParameters = Acts::CurvilinearTrackParameters;

///// Container of reconstructed track states for multiple tracks.
// using TrackParametersContainer = std::vector<TrackParameters>;

///// MultiTrajectory definition
// using Trajectory = Acts::MultiTrajectory<SourceLink>;

///// Container for the truth fitting/finding track(s)
// using TrajectoryContainer = std::vector<SimMultiTrajectory>;

namespace Jug::Reco {

/** Initial Track parameters from MC truth.
 *
 *  TrackParmetersContainer
 *
 *  \ingroup tracking
 */
class TrackParamClusterInit : public GaudiAlgorithm {
private:
  using Clusters = edm4eic::ClusterCollection;

  DataHandle<Clusters> m_inputClusters{"inputClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<TrackParametersContainer> m_outputInitialTrackParameters{"outputInitialTrackParameters",
                                                                      Gaudi::DataHandle::Writer, this};

public:
  TrackParamClusterInit(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputClusters", m_inputClusters, "Input clusters");
    declareProperty("outputInitialTrackParameters", m_outputInitialTrackParameters, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
    if (randSvc == nullptr) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collection
    const auto* const clusters = m_inputClusters.get();
    // Create output collections
    auto* init_trk_params = m_outputInitialTrackParameters.createAndPut();

    for (const auto& c : *clusters) {

      using Acts::UnitConstants::GeV;
      using Acts::UnitConstants::MeV;
      using Acts::UnitConstants::mm;
      using Acts::UnitConstants::ns;

      double p = c.getEnergy() * GeV;
      if (p < 0.1 * GeV) {
        continue;
      }
      double len    = edm4eic::magnitude(c.getPosition());
      auto momentum = c.getPosition() * p / len;

      Acts::BoundVector params;
      params(Acts::eBoundLoc0)   = 0.0 * mm;
      params(Acts::eBoundLoc1)   = 0.0 * mm;
      params(Acts::eBoundPhi)    = edm4eic::angleAzimuthal(momentum);
      params(Acts::eBoundTheta)  = edm4eic::anglePolar(momentum);
      params(Acts::eBoundQOverP) = 1 / p;
      params(Acts::eBoundTime)   = 0 * ns;

      auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0, 0, 0});

      if (msgLevel(MSG::DEBUG)) {
        debug() << "Invoke track finding seeded by truth particle with p = " << p / GeV << " GeV" << endmsg;
      }

      // add both charges to the track candidate...
      init_trk_params->push_back({pSurface, params, 1});

      Acts::BoundVector params2;
      params2(Acts::eBoundLoc0)   = 0.0 * mm;
      params2(Acts::eBoundLoc1)   = 0.0 * mm;
      params2(Acts::eBoundPhi)    = edm4eic::angleAzimuthal(momentum);
      params2(Acts::eBoundTheta)  = edm4eic::anglePolar(momentum);
      params2(Acts::eBoundQOverP) = -1 / p;
      params2(Acts::eBoundTime)   = 0 * ns;
      init_trk_params->push_back({pSurface, params2, -1});

      // acts v1.2.0:
      // init_trk_params->emplace_back(Acts::Vector4(0 * mm, 0 * mm, 0 * mm, 0),
      //                              Acts::Vector3(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, -1,
      //                              std::make_optional(cov));
      // debug() << init_trk_params->back() << endmsg;
      // init_trk_params->emplace_back(Acts::Vector4(0 * mm, 0 * mm, 0 * mm, 0),
      //                              Acts::Vector3(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, 1,
      //                              std::make_optional(cov));
      ////debug() << init_trk_params->back() << endmsg;
      // init_trk_params->emplace_back(Acts::Vector4(0 * mm, 0 * mm, 0 * mm, 0),
      //                              Acts::Vector3(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, 0,
      //                              std::make_optional(cov));
    }
    return StatusCode::SUCCESS;
  }
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(TrackParamClusterInit)

} // namespace Jug::Reco
