// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wenqing Fan, Barak Schmookler, Whitney Armstrong, Sylvester Joosten

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
#include "edm4eic/TrackerHitCollection.h"
#include "edm4eic/TrackParametersCollection.h"
#include "edm4eic/TrajectoryCollection.h"
#include "edm4eic/TrackSegmentCollection.h"
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
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"

#include "edm4eic/vector_utils.h"

#include <cmath>

using BoundTrackParamPtr = std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;


namespace Jug::Reco {

  /** Extrac the particles form fit trajectories.
   *
   * \ingroup tracking
   */
   class TrackProjector : public GaudiAlgorithm {
   private:
    DataHandle<TrajectoriesContainer>        m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<edm4eic::TrackSegmentCollection> m_outputTrackSegments{"outputTrackSegments", Gaudi::DataHandle::Writer, this};

    Gaudi::Property<unsigned int> m_firstInVolumeID{this, "firstInVolumeID", 0};
    Gaudi::Property<std::string> m_firstInVolumeName{this, "firstInVolumeName", ""};
    Gaudi::Property<float> m_firstSmallerThanZ{this, "firstSmallerThanZ", 0};
    Gaudi::Property<float> m_firstGreaterThanZ{this, "firstGreaterThanZ", 0};
    Gaudi::Property<float> m_firstGreaterThanR{this, "firstGreaterThanR", -1};

    SmartIF<IGeoSvc> m_geoSvc;
    Acts::GeometryContext m_geoContext;
    Acts::MagneticFieldContext m_fieldContext;

    std::shared_ptr<Acts::DiscSurface> hcalEndcapNSurf;

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

      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }

      auto transform = Acts::Transform3::Identity();
	
      // make a reference disk to mimic electron-endcap HCal
      const auto hcalEndcapNZ = -3322.;
      const auto hcalEndcapNMinR = 83.01;
      const auto hcalEndcapNMaxR = 950;
      auto hcalEndcapNBounds = std::make_shared<Acts::RadialBounds>(hcalEndcapNMinR, hcalEndcapNMaxR);
      auto hcalEndcapNTrf = transform * Acts::Translation3(Acts::Vector3(0, 0, hcalEndcapNZ));
      hcalEndcapNSurf = Acts::Surface::makeShared<Acts::DiscSurface>(hcalEndcapNTrf, hcalEndcapNBounds);

      return StatusCode::SUCCESS;
    }

    BoundTrackParamPtrResult propagateTrack(const Acts::BoundTrackParameters& params, const SurfacePtr& targetSurf) {

	    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = m_geoSvc->trackingGeometry();
	    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField = m_geoSvc->getFieldProvider();
	    using Stepper            = Acts::EigenStepper<>;
	    using Propagator         = Acts::Propagator<Stepper>;
	    Stepper stepper(magneticField);
	    Propagator propagator(stepper);
	    // Acts::Logging::Level logLevel = Acts::Logging::FATAL
	    Acts::Logging::Level logLevel = Acts::Logging::DEBUG;
	    
	    ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("ProjectTrack Logger", logLevel));
      
	    Acts::PropagatorOptions<> options(m_geoContext,m_fieldContext,Acts::LoggerWrapper{logger()});

	    auto result = propagator.propagate(params, *targetSurf, options);
   
	    if(result.ok()) return std::move((*result).endParameters);
	    return result.error();
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

        edm4eic::MutableTrackSegment track_segment;

        // visit the track points
        mj.visitBackwards(trackTip, [&](auto&& trackstate) {
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
          const decltype(edm4eic::TrackPoint::position) position {
            static_cast<float>(global.x()),
            static_cast<float>(global.y()),
            static_cast<float>(global.z())
          };

          // local position
          const decltype(edm4eic::TrackParametersData::loc) loc {
            static_cast<float>(parameter[Acts::eBoundLoc0]),
            static_cast<float>(parameter[Acts::eBoundLoc1])
          };
          const decltype(edm4eic::TrackParametersData::locError) locError {
            static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)),
            static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)),
            static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc1))
          };
          const decltype(edm4eic::TrackPoint::positionError) positionError{0, 0, 0};
          const decltype(edm4eic::TrackPoint::momentum) momentum = edm4eic::sphericalToVector(
            static_cast<float>(1.0 / std::abs(parameter[Acts::eBoundQOverP])),
            static_cast<float>(parameter[Acts::eBoundTheta]),
            static_cast<float>(parameter[Acts::eBoundPhi])
          );
          const decltype(edm4eic::TrackPoint::momentumError) momentumError {
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
          const decltype(edm4eic::TrackPoint::directionError) directionError {
            static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
            static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundPhi)),
            static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundPhi))
          };
          const float pathLength = static_cast<float>(trackstate.pathLength());
          const float pathLengthError = 0;

          // Store track point
          track_segment.addToPoints({
            position,
            positionError,
            momentum,
            momentumError,
            time,
            timeError,
            theta,
            phi,
            directionError,
            pathLength,
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

	      //=================================================
	      //Track projection
	      //Reference sPHENIX code: https://github.com/sPHENIX-Collaboration/coresoftware/blob/335e6da4ccacc8374cada993485fe81d82e74a4f/offline/packages/trackreco/PHActsTrackProjection.h
	      //=================================================
	      const auto& boundParam = traj.trackParameters(trackTip);	
	
	      // project track parameters to electron endcap hcal surface
	      auto result = propagateTrack(boundParam, hcalEndcapNSurf);
              if(result.ok()){
		      auto trackStateParams = std::move(**result);
                      auto projectionPos = trackStateParams.position(m_geoContext);

		      if (msgLevel(MSG::DEBUG)) {
         		 debug() << "X projection is " << projectionPos(0) << endmsg;
			 debug() << "Y projection is " << projectionPos(1) << endmsg;
			 debug() << "Z projection is " << projectionPos(2) << endmsg;
        	      }

		      const decltype(edm4eic::TrackPoint::position) proj_position {
		           static_cast<float>(projectionPos(0)),
            		   static_cast<float>(projectionPos(1)),
            		   static_cast<float>(projectionPos(2))
          	      };
		      const decltype(edm4eic::TrackPoint::positionError) proj_positionError{0, 0, 0};
		      const decltype(edm4eic::TrackPoint::momentum) proj_momentum{0, 0, 0};
		      const decltype(edm4eic::TrackPoint::momentumError) proj_momentumError{0, 0, 0};
		      const float proj_time = 0;
		      const float proj_timeError = 0;
		      const float proj_theta = 0;
                      const float proj_phi = 0;
		      const decltype(edm4eic::TrackPoint::directionError) proj_directionError{0, 0, 0};
		      const float proj_pathLength = 0;
		      const float proj_pathLengthError = 0;

		      // Store projection point
		      track_segment.addToPoints({
            		proj_position,
          		proj_positionError,
            		proj_momentum,
            		proj_momentumError,
            		proj_time,
            		proj_timeError,
            		proj_theta,
           	 	proj_phi,
            		proj_directionError,
            		proj_pathLength,
            		proj_pathLengthError
          	     });

	      }

	// Set associated track
	// track_segment.setTrack(traj);	

        // Add to output collection
        track_segments->push_back(track_segment);
      }

      return StatusCode::SUCCESS;
    }

  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(TrackProjector)

} // namespace Jug::Reco
