// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

#ifndef JUGGLER_JUGRECO_CKFTracking_HH
#define JUGGLER_JUGRECO_CKFTracking_HH

#include <functional>
#include <random>
#include <stdexcept>
#include <vector>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>
#include "JugTrack/IActsGeoSvc.h"
#include "JugTrack/DD4hepBField.h"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"

#include "edm4eic/TrackerHitCollection.h"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"

namespace Jug::Reco {

/** Fitting algorithm implmentation .
 *
 * \ingroup tracking
 */
class CKFTracking : public GaudiAlgorithm {
public:
  /// Track finder function that takes input measurements, initial trackstate
  /// and track finder options and returns some track-finder-specific result.
#if Acts_VERSION_MAJOR >= 36
  using TrackFinderOptions =
      Acts::CombinatorialKalmanFilterOptions<ActsExamples::IndexSourceLinkAccessor::Iterator,
                                             ActsExamples::TrackContainer>;
#else
  using TrackFinderOptions =
      Acts::CombinatorialKalmanFilterOptions<ActsExamples::IndexSourceLinkAccessor::Iterator,
                                             Acts::VectorMultiTrajectory>;
#endif
  using TrackFinderResult =
      Acts::Result<std::vector<ActsExamples::TrackContainer::TrackProxy>>;

  /// Find function that takes the above parameters
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class CKFTrackingFunction {
   public:
    virtual ~CKFTrackingFunction() = default;
    virtual TrackFinderResult operator()(const ActsExamples::TrackParameters&,
                                         const TrackFinderOptions&,
                                         ActsExamples::TrackContainer&) const = 0;
  };

  /// Create the track finder function implementation.
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static std::shared_ptr<CKFTrackingFunction> makeCKFTrackingFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField);

public:
#if Acts_VERSION_MAJOR < 37 || (Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR < 1)
  DataHandle<ActsExamples::IndexSourceLinkContainer> m_inputSourceLinks{"inputSourceLinks", Gaudi::DataHandle::Reader, this};
#endif
  DataHandle<ActsExamples::MeasurementContainer> m_inputMeasurements{"inputMeasurements", Gaudi::DataHandle::Reader, this};
  DataHandle<ActsExamples::TrackParametersContainer> m_inputInitialTrackParameters{"inputInitialTrackParameters",
                                                                     Gaudi::DataHandle::Reader, this};
  DataHandle<ActsExamples::ConstTrackContainer> m_outputTracks{"outputTracks", Gaudi::DataHandle::Writer, this};
  DataHandle<ActsExamples::TrajectoriesContainer> m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<std::vector<double>> m_etaBins{this, "etaBins", {}};
  Gaudi::Property<std::vector<double>> m_chi2CutOff{this, "chi2CutOff", {15.}};
  Gaudi::Property<std::vector<size_t>> m_numMeasurementsCutOff{this, "numMeasurementsCutOff", {10}};

  std::shared_ptr<CKFTrackingFunction> m_trackFinderFunc;
  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IActsGeoSvc> m_actsGeoSvc;

  std::shared_ptr<const Jug::BField::DD4hepBField> m_BField = nullptr;
  Acts::GeometryContext m_geoctx;
  Acts::CalibrationContext m_calibctx;
  Acts::MagneticFieldContext m_fieldctx;

  Acts::MeasurementSelector::Config m_sourcelinkSelectorCfg;
  Acts::Logging::Level m_actsLoggingLevel = Acts::Logging::INFO;

  CKFTracking(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  StatusCode execute() override;
};

} // namespace Jug::Reco

#endif
