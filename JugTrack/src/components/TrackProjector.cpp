// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 wfan, Whitney Armstrong, Sylvester Joosten

#include <algorithm>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

// Event Model related classes
#include "eicd/BasicParticleCollection.h"
#include "eicd/TrackerHitCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/TrajectoryCollection.h"
#include "eicd/TrackSegmentCollection.h"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Trajectories.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

#include "eicd/vector_utils.h"

#include <cmath>

namespace Jug::Reco {

  /** Extrac the particles form fit trajectories.
   *
   * \ingroup tracking
   */
   class TrackProjector : public GaudiAlgorithm {
   private:
    DataHandle<TrajectoriesContainer>        m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<eicd::TrackSegmentCollection> m_outputTrackSegments{"outputTrackSegments", Gaudi::DataHandle::Writer, this};

    Gaudi::Property<unsigned int> m_firstInVolumeID{this, "firstInVolumeID", 0};
    Gaudi::Property<std::string> m_firstInVolumeName{this, "firstInVolumeName", ""};
    Gaudi::Property<float> m_firstSmallerThanZ{this, "firstSmallerThanZ", 0};
    Gaudi::Property<float> m_firstGreaterThanZ{this, "firstGreaterThanZ", 0};
    Gaudi::Property<float> m_firstGreaterThanR{this, "firstGreaterThanR", -1};

    Acts::GeometryContext m_geoContext;

    public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    TrackProjector(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputTrajectories", m_inputTrajectories,"");
          declareProperty("outputTrackSegments", m_outputTrackSegments, "");
        }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const auto* const trajectories = m_inputTrajectories.get();
      // create output collections
      auto* track_segments = m_outputTrackSegments.createAndPut();

      if (msgLevel(MSG::DEBUG)) {
        debug() << std::size(*trajectories) << " trajectories " << endmsg;
      }

      // Loop over the trajectories
      for (const auto& traj : *trajectories) {
        // Get the entry index for the single trajectory
        // The trajectory entry indices and the multiTrajectory
        const auto& mj        = traj.multiTrajectory();
        const auto& trackTips = traj.tips();
        if (msgLevel(MSG::DEBUG)) {
          debug() << "# of elements in trackTips " <<trackTips.size() << endmsg;
        }

        // Skip empty
        if (trackTips.empty()) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "Empty multiTrajectory." << endmsg;
          }
          continue;
        }
        auto& trackTip = trackTips.front();

        // Collect the trajectory summary info
        auto trajState       = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
        int  m_nMeasurements = trajState.nMeasurements;
        int  m_nStates       = trajState.nStates;
        int  m_nCalibrated   = 0;
        if (msgLevel(MSG::DEBUG)) {
          debug() << "n measurement in trajectory " << m_nMeasurements << endmsg;
          debug() << "n state in trajectory " << m_nStates << endmsg;
        }

        eicd::MutableTrackSegment track_segment;

        // visit the track points
        mj.visitBackwards(trackTip, [&](auto&& trackstate) {
          auto pathLength = trackstate.pathLength();

          // get volume info
          auto geoID = trackstate.referenceSurface().geometryId();
          auto volume = geoID.volume();
          auto layer = geoID.layer();
          if (trackstate.hasCalibrated()) {
            m_nCalibrated++;
          }

          // get track state parameters and their covariances
          const auto& parameter = trackstate.predicted();
          const auto& covariance = trackstate.predictedCovariance();

          // convert local to global
          auto global = trackstate.referenceSurface().localToGlobal(
            m_geoContext,
            {parameter[Acts::eBoundLoc0], parameter[Acts::eBoundLoc1]},
            {0, 0, 0}
          );
          // global position
          const eicd::Vector3f position {
            global.x(),
            global.y(),
            global.z()
          };

          // local position
          const eicd::Vector2f loc {
            parameter[Acts::eBoundLoc0],
            parameter[Acts::eBoundLoc1]
          };
          const eicd::Cov2f covLoc {
            static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)),
            static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)),
            static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc1))
          };
          const eicd::Cov3f positionError{0, 0, 0};
          const eicd::Vector3f momentum = eicd::sphericalToVector(
            1.0 / std::abs(parameter[Acts::eBoundQOverP]),
            parameter[Acts::eBoundTheta],
            parameter[Acts::eBoundPhi]
          );
          const eicd::Cov3f covMomentum {
            static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
            static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundPhi)),
            static_cast<float>(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)),
            static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundPhi)),
            static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundQOverP)),
            static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundQOverP))
          };
          const float time{static_cast<float>(parameter(Acts::eBoundTime))};
          const float timeError{sqrt(static_cast<float>(covariance(Acts::eBoundTime, Acts::eBoundTime)))};
          const float theta(parameter[Acts::eBoundTheta]);
          const float phi(parameter[Acts::eBoundPhi]);
          const eicd::Cov2f directionError {
            static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
            static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundPhi)),
            static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundPhi))
          };
          const float pathLengthError = 0;

          // Store track point
          track_segment.addToPoints({
            position,
            positionError,
            momentum,
            covMomentum,
            time,
            timeError,
            theta,
            phi,
            directionError,
            static_cast<float>(pathLength),
            pathLengthError
          });

          if (msgLevel(MSG::DEBUG)) {
            debug() << "******************************" << endmsg;
            debug() << "predicted variables: \n" << trackstate.predicted() << endmsg;
            debug() << "geoID = " << geoID << endmsg;
            debug() << "volume = " << volume << ", layer = " << layer << endmsg;
            debug() << "pathlength = " << pathLength << endmsg;
            debug() << "hasCalibrated = " << trackstate.hasCalibrated() << endmsg;
            debug() << "******************************" << endmsg;
          }
        });

        if (msgLevel(MSG::DEBUG)) {
          debug() << "n calibrated state in trajectory " << m_nCalibrated << endmsg;
        }

        // Add to output collection
        track_segments->push_back(track_segment);
      }

      return StatusCode::SUCCESS;
    }

  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(TrackProjector)

} // namespace Jug::Reco
