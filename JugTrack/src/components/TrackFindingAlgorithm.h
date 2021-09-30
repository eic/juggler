#ifndef JUGGLER_JUGRECO_TrackFindingAlgorithm_HH
#define JUGGLER_JUGRECO_TrackFindingAlgorithm_HH

#include <functional>
#include <random>
#include <stdexcept>
#include <vector>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/BField/DD4hepBField.h"
#include "JugTrack/GeometryContainers.hpp"
#include "JugTrack/Index.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Trajectories.hpp"

#include "eicd/TrackerHitCollection.h"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"

namespace Jug::Reco {

/** Fitting algorithm implmentation .
 *
 * \ingroup tracking
 */
class TrackFindingAlgorithm : public GaudiAlgorithm {
public:
  /// Track finder function that takes input measurements, initial trackstate
  /// and track finder options and returns some track-finder-specific result.
  using TrackFinderOptions  = Acts::CombinatorialKalmanFilterOptions<IndexSourceLinkAccessor, MeasurementCalibrator, Acts::MeasurementSelector>;
  using TrackFinderResult   = std::vector<Acts::Result<Acts::CombinatorialKalmanFilterResult<IndexSourceLink>>>;
  using TrackFinderFunction = std::function<TrackFinderResult(
      const IndexSourceLinkContainer&, const TrackParametersContainer&, const TrackFinderOptions&)>;

  /// Create the track finder function implementation.
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static TrackFinderFunction makeTrackFinderFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                                     std::shared_ptr<const Acts::MagneticFieldProvider> magneticField);

public:
  DataHandle<IndexSourceLinkContainer> m_inputSourceLinks{"inputSourceLinks", Gaudi::DataHandle::Reader, this};
  DataHandle<MeasurementContainer> m_inputMeasurements{"inputMeasurements", Gaudi::DataHandle::Reader, this};
  DataHandle<TrackParametersContainer> m_inputInitialTrackParameters{"inputInitialTrackParameters",
                                                                     Gaudi::DataHandle::Reader, this};
  DataHandle<TrajectoriesContainer> m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<double> m_chi2CutOff{this, "chi2CutOff", 15.};
  Gaudi::Property<size_t> m_numMeasurementsCutOff{this, "numMeasurementsCutOff", 10};

  TrackFinderFunction m_trackFinderFunc;
  SmartIF<IGeoSvc> m_geoSvc;

  std::shared_ptr<const Jug::BField::DD4hepBField> m_BField = nullptr;
  Acts::GeometryContext m_geoctx;
  Acts::CalibrationContext m_calibctx;
  Acts::MagneticFieldContext m_fieldctx;

  Acts::MeasurementSelector::Config m_sourcelinkSelectorCfg;
  Acts::Logging::Level m_actsLoggingLevel = Acts::Logging::INFO;

  TrackFindingAlgorithm(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  StatusCode execute() override;
};

} // namespace Jug::Reco

#endif
