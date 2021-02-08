#ifndef JUGGLER_JUGRECO_TrackFindingAlgorithm_HH
#define JUGGLER_JUGRECO_TrackFindingAlgorithm_HH

#include "JugReco/GeometryContainers.hpp"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"


#include <functional>
#include <random>
#include <stdexcept>
#include <vector>

#include "JugReco/Index.hpp"
#include "JugReco/IndexSourceLink.hpp"
#include "JugReco/Track.hpp"
#include "JugReco/BField.h"
#include "JugReco/Measurement.hpp"

#include "eicd/TrackerHitCollection.h"

//#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/Definitions/Common.hpp"


namespace Jug::Reco {

  class TrackFindingAlgorithm : public GaudiAlgorithm {
  public:
    /// Track finder function that takes input measurements, initial trackstate
    /// and track finder options and returns some track-finder-specific result.
    using TrackFinderOptions = Acts::CombinatorialKalmanFilterOptions<MeasurementCalibrator, Acts::MeasurementSelector>;
    using TrackFinderResult  = std::vector<Acts::Result<Acts::CombinatorialKalmanFilterResult<IndexSourceLink>>>;
    using TrackFinderFunction = std::function<TrackFinderResult(
        const IndexSourceLinkContainer&, const TrackParametersContainer&, const TrackFinderOptions&)>;

    /// Create the track finder function implementation.
    /// The magnetic field is intentionally given by-value since the variant
    /// contains shared_ptr anyways.
    static TrackFinderFunction makeTrackFinderFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                                       BFieldVariant                        magneticField);

  public:
    DataHandle<IndexSourceLinkContainer> m_inputSourceLinks{"inputSourceLinks", Gaudi::DataHandle::Reader, this};
    DataHandle<MeasurementContainer>     m_inputMeasurements{"inputMeasurements", Gaudi::DataHandle::Reader, this};
    DataHandle<TrackParametersContainer> m_inputInitialTrackParameters{"inputInitialTrackParameters",
                                                                       Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoryContainer>      m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};
    TrackFinderFunction                  m_trackFinderFunc;
    SmartIF<IGeoSvc>                     m_geoSvc;

    std::shared_ptr<Acts::ConstantBField> m_BField = nullptr;
    Acts::GeometryContext                 m_geoctx;
    Acts::CalibrationContext              m_calibctx;
    Acts::MagneticFieldContext            m_fieldctx;

    Acts::MeasurementSelector::Config m_sourcelinkSelectorCfg;

    TrackFindingAlgorithm(const std::string& name, ISvcLocator* svcLoc);

    StatusCode initialize() override;

    StatusCode execute() override;
  };

} // namespace Jug::Reco

#endif
