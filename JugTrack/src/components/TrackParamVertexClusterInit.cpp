// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten

#include <cmath>

// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugTrack/Track.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "eicd/ClusterCollection.h"
#include "eicd/TrackerHitCollection.h"
#include "eicd/vector_utils.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

/** Initial Track parameters from clusters and vertex hits.
 *
 *
 * The momentum of the initial track is estimated from the cluster  energy and
 * the direction is set using the vertex hits.
 *
 * \ingroup tracking
 */
class TrackParamVertexClusterInit : public GaudiAlgorithm {
private:
  using Clusters   = eicd::ClusterCollection;
  using VertexHits = eicd::TrackerHitCollection;

  DataHandle<VertexHits> m_inputVertexHits{"inputVertexHits", Gaudi::DataHandle::Reader, this};
  DataHandle<Clusters> m_inputClusters{"inputClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<TrackParametersContainer> m_outputInitialTrackParameters{"outputInitialTrackParameters",
                                                                      Gaudi::DataHandle::Writer, this};
  Gaudi::Property<double> m_maxHitRadius{this, "maxHitRadius", 40.0 * mm};

public:
  TrackParamVertexClusterInit(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputVertexHits", m_inputVertexHits, "Vertex tracker hits");
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
    const auto* const vtx_hits = m_inputVertexHits.get();
    // Create output collections
    auto* init_trk_params = m_outputInitialTrackParameters.createAndPut();

    double max_radius = m_maxHitRadius.value();

    for (const auto& c : *clusters) {

      using Acts::UnitConstants::GeV;
      using Acts::UnitConstants::MeV;
      using Acts::UnitConstants::mm;
      using Acts::UnitConstants::ns;

      double p_cluster = c.getEnergy() * GeV;
      if (p_cluster / GeV < 0.1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << " skipping cluster with energy " << p_cluster / GeV << " GeV" << endmsg;
        }
        continue;
      }

      auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0, 0, 0});
      for (const auto& t : *vtx_hits) {

        double len = std::hypot(t.getPosition().x, t.getPosition().y, t.getPosition().z);
        if (len > max_radius) {
          continue;
        }

        auto momentum = t.getPosition() * p_cluster / len;

        Acts::BoundVector params;
        params(Acts::eBoundLoc0)   = 0.0 * mm;
        params(Acts::eBoundLoc1)   = 0.0 * mm;
        params(Acts::eBoundPhi)    = eicd::angleAzimuthal(momentum);
        params(Acts::eBoundTheta)  = eicd::anglePolar(momentum);
        params(Acts::eBoundQOverP) = 1 / p_cluster;
        params(Acts::eBoundTime)   = 0 * ns;

        // debug() << "Invoke track finding seeded by truth particle with p = " << p / GeV << " GeV" << endmsg;

        // add both charges to the track candidate...
        init_trk_params->push_back({pSurface, params, 1});

        Acts::BoundVector params2;
        params2(Acts::eBoundLoc0)   = 0.0 * mm;
        params2(Acts::eBoundLoc1)   = 0.0 * mm;
        params2(Acts::eBoundPhi)    = eicd::angleAzimuthal(momentum);
        params2(Acts::eBoundTheta)  = eicd::anglePolar(momentum);
        params2(Acts::eBoundQOverP) = -1 / p_cluster;
        params2(Acts::eBoundTime)   = 0 * ns;
        init_trk_params->push_back({pSurface, params2, -1});
      }
      // init_trk_params->emplace_back(Acts::Vector4(0 * mm, 0 * mm, 0 * mm, 0),
      //                              Acts::Vector3(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, 1,
      //                              std::make_optional(cov));
      // init_trk_params->emplace_back(Acts::Vector4(0 * mm, 0 * mm, 0 * mm, 0),
      //                              Acts::Vector3(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, 0,
      //                              std::make_optional(cov));
      // debug() << "Invoke track finding seeded by truth particle with p = " << p_cluster/GeV  << " GeV" << endmsg;
    }
    return StatusCode::SUCCESS;
  }
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(TrackParamVertexClusterInit)

} // namespace Jug::Reco
