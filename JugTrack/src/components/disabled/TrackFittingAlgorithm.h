// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#ifndef JUGGLER_JUGRECO_TrackFittingAlgorithm_HH
#define JUGGLER_JUGRECO_TrackFittingAlgorithm_HH 1

#include <functional>
#include <stdexcept>
#include <vector>
#include <random>
#include <stdexcept>

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "Gaudi/Property.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>
#include <k4ActsTracking/IActsGeoSvc.h>
#include "JugBase/BField/DD4hepBField.h"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"

#include "edm4eic/TrackerHitCollection.h"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Definitions/Common.hpp"


namespace Jug::Reco {

  /** Fitting algorithm implmentation .
   *
   * \ingroup tracking
   */
  class TrackFittingAlgorithm : public Gaudi::Algorithm {
  public:
    /// Track fitter function that takes input measurements, initial trackstate
    /// and fitter options and returns some track-fitter-specific result.
    using TrackFitterOptions =
        Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory>;

    using FitterResult =
	Acts::Result<Acts::KalmanFitterResult<Acts::VectorMultiTrajectory>>;

    /// Fit function that takes input measurements, initial trackstate and fitter
    using FitterFunction = std::function<FitterResult(
      const std::vector<ActsExamples::IndexSourceLink>&, const ActsExamples::TrackParameters&, const TrackFitterOptions&)>;

  public:
    mutable DataHandle<const ActsExamples::IndexSourceLinkContainer> m_inputSourceLinks{"inputSourceLinks", Gaudi::DataHandle::Reader, this};
    DataHandle<ActsExamples::TrackParametersContainer> m_initialTrackParameters{"initialTrackParameters", Gaudi::DataHandle::Reader, this};
    mutable DataHandle<const ActsExamples::MeasurementContainer>     m_inputMeasurements{"inputMeasurements", Gaudi::DataHandle::Reader, this};
    mutable DataHandle<const ActsExamples::ProtoTrackContainer>      m_inputProtoTracks{"inputProtoTracks", Gaudi::DataHandle::Reader, this};
    DataHandle<ActsExamples::TrajectoriesContainer>    m_foundTracks{"foundTracks", Gaudi::DataHandle::Reader, this};
    mutable DataHandle<ActsExamples::TrajectoriesContainer>    m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};

    FitterFunction                        m_trackFittingFunc;
    SmartIF<IGeoSvc>                      m_geoSvc;
    SmartIF<IActsGeoSvc>                  m_actsGeoSvc;
    std::shared_ptr<const Jug::BField::DD4hepBField> m_BField = nullptr;
    Acts::GeometryContext                 m_geoctx;
    Acts::CalibrationContext              m_calibctx;
    Acts::MagneticFieldContext            m_fieldctx;

    //Acts::CKFSourceLinkSelector::Config m_sourcelinkSelectorCfg;

    TrackFittingAlgorithm(const std::string& name, ISvcLocator* svcLoc);

    /** Create the track finder function implementation.
     *  The magnetic field is intentionally given by-value since the variant
     *  contains shared_ptr anyways.
     */
    static FitterFunction
    makeTrackFittingFunction(std::shared_ptr<const Acts::TrackingGeometry>      trackingGeometry,
                             std::shared_ptr<const Acts::MagneticFieldProvider> magneticField);

    StatusCode initialize() override;

    StatusCode execute(const EventContext&) const override;
   private:
    /// Helper function to call correct FitterFunction
    FitterResult fitTrack(
        const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
        const ActsExamples::TrackParameters& initialParameters,
        const TrackFitterOptions& options
        ) const;
  };

  inline TrackFittingAlgorithm::FitterResult
  TrackFittingAlgorithm::fitTrack(const std::vector<ActsExamples::IndexSourceLink>& sourceLinks,
                                  const ActsExamples::TrackParameters&              initialParameters,
                                  const TrackFitterOptions&           options) const
  {
    return m_trackFittingFunc(sourceLinks, initialParameters, options);
  }

} // namespace Jug::Reco

#endif
