// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten

#include <cmath>
// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>
#include "ActsExamples/EventData/Track.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"

#include "edm4eic/TrackerHitCollection.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4hep/utils/vector_utils.h"

#include "Acts/Surfaces/PerigeeSurface.hpp"

  ///// (Reconstructed) track parameters e.g. close to the vertex.
  //using TrackParameters = Acts::CurvilinearTrackParameters;

  ///// Container of reconstructed track states for multiple tracks.
  //using TrackParametersContainer = std::vector<TrackParameters>;

  ///// MultiTrajectory definition
  //using Trajectory = Acts::MultiTrajectory<SourceLink>;

  ///// Container for the truth fitting/finding track(s)
  //using TrajectoryContainer = std::vector<SimMultiTrajectory>;

namespace Jug::Reco {

  /** Initial Track parameters from MC truth.
   *
   *  TrackParmetersContainer
   *
   *  \ingroup tracking
   */
  class TrackParamImagingClusterInit : public Gaudi::Algorithm {
  private:
    mutable DataHandle<edm4eic::ClusterCollection>          m_inputClusters{"inputClusters", Gaudi::DataHandle::Reader, this};
    mutable DataHandle<ActsExamples::TrackParametersContainer> m_outputInitialTrackParameters{"outputInitialTrackParameters",
                                                                        Gaudi::DataHandle::Writer, this};

  public:
    TrackParamImagingClusterInit(const std::string& name, ISvcLocator* svcLoc)
        : Gaudi::Algorithm(name, svcLoc) {
      declareProperty("inputClusters", m_inputClusters, "Input clusters");
      declareProperty("outputInitialTrackParameters", m_outputInitialTrackParameters, "");
    }

    StatusCode initialize() override {
      if (Gaudi::Algorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      IRndmGenSvc* randSvc = Gaudi::svcLocator()->service<IRndmGenSvc>("RndmGenSvc", true);
      if (randSvc == nullptr) {
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute(const EventContext&) const override {
      // input collection
      const auto* const clusters = m_inputClusters.get();
      // Create output collections
      auto* init_trk_params = m_outputInitialTrackParameters.createAndPut();

      for(const auto& c : *clusters) {

        using Acts::UnitConstants::GeV;
        using Acts::UnitConstants::MeV;
        using Acts::UnitConstants::mm;
        using Acts::UnitConstants::ns;

        const double p = c.getEnergy()*GeV;
        // FIXME hardcoded value
        if( p < 0.1*GeV) {
          continue;
        }
        const double theta = edm4hep::utils::anglePolar(c.getPosition());
        const double phi = edm4hep::utils::angleAzimuthal(c.getPosition());

        Acts::BoundVector  params;
        params(Acts::eBoundLoc0)   = 0.0 * mm ;
        params(Acts::eBoundLoc1)   = 0.0 * mm ;
        params(Acts::eBoundPhi)    = phi;
        params(Acts::eBoundTheta)  = theta;
        params(Acts::eBoundQOverP) = 1/p;
        params(Acts::eBoundTime)   = 0 * ns;

        auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0,0,0});

        debug() << "Invoke track finding seeded by truth particle with p = " << p/GeV  << " GeV" << endmsg;

        // add both charges to the track candidate...
        init_trk_params->push_back({pSurface, params, {}, Acts::ParticleHypothesis::pion()});

        Acts::BoundVector  params2;
        params2(Acts::eBoundLoc0)   = 0.0 * mm ;
        params2(Acts::eBoundLoc1)   = 0.0 * mm ;
        params2(Acts::eBoundPhi)    = phi;
        params2(Acts::eBoundTheta)  = theta;
        params2(Acts::eBoundQOverP) = -1/p;
        params2(Acts::eBoundTime)   = 0 * ns;
        init_trk_params->push_back({pSurface, params2, {}, Acts::ParticleHypothesis::pion()});

      }
      return StatusCode::SUCCESS;
    }
  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(TrackParamImagingClusterInit)

} // namespace Jug::reco
