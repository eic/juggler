#ifndef JUGGLER_JUGRECO_TrackFindingAlgorithm_HH
#define JUGGLER_JUGRECO_TrackFindingAlgorithm_HH

#include "JugReco/GeometryContainers.hpp"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
//#include "GaudiAlg/Transformer.h"
//#include "GaudiAlg/GaudiTool.h"
//#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

//#include "Acts/Geometry/TrackingGeometry.hpp"
//#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
//#include "Acts/Utilities/Definitions.hpp"
//#include "Acts/Utilities/Helpers.hpp"
//#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <stdexcept>
#include <vector>

#include "JugReco/SourceLinks.h"
#include "JugReco/Track.hpp"
#include "JugReco/BField.h"

#include "eicd/TrackerHitCollection.h"

//#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFinder/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinder/CombinatorialKalmanFilter.hpp"

//#include "Acts/Fitter/GainMatrixSmoother.hpp"
//#include "Acts/Fitter/GainMatrixUpdater.hpp"
//#include "Acts/MagneticField/ConstantBField.hpp"
//#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
//#include "Acts/MagneticField/SharedBField.hpp"
//#include "Acts/Propagator/EigenStepper.hpp"
//#include "Acts/Propagator/Navigator.hpp"
//#include "Acts/Propagator/Propagator.hpp"
//#include "Acts/Utilities/Units.hpp"

#include <random>
#include <stdexcept>

namespace Jug::Reco {

  class TrackFindingAlgorithm : public GaudiAlgorithm {
  public:
    using TrackFinderResult = Acts::Result<Acts::CombinatorialKalmanFilterResult<SourceLink>>;

    /// Track finding function that takes input measurements, initial trackstate
    /// and track finder options and returns some track-finding-specific result.
    using CKFOptions = Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>;

    using TrackFinderFunction =
        std::function<TrackFinderResult(const SourceLinkContainer&, const TrackParameters&, const CKFOptions&)>;

  public:
    DataHandle<SourceLinkContainer>      m_inputSourceLinks{"inputSourceLinks", Gaudi::DataHandle::Reader, this};
    DataHandle<TrackParametersContainer> m_inputInitialTrackParameters{"inputInitialTrackParameters",
                                                                       Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoryContainer>      m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};
    TrackFinderFunction                  m_trackFinderFunc;
    SmartIF<IGeoSvc>                     m_geoSvc;

    std::shared_ptr<Acts::ConstantBField> m_BField = nullptr;
    Acts::GeometryContext                 m_geoctx;
    Acts::CalibrationContext              m_calibctx;
    Acts::MagneticFieldContext            m_fieldctx;

    Acts::CKFSourceLinkSelector::Config m_sourcelinkSelectorCfg;

    TrackFindingAlgorithm(const std::string& name, ISvcLocator* svcLoc);

    /** Create the track finder function implementation.
     *  The magnetic field is intentionally given by-value since the variant
     *  contains shared_ptr anyways.
     */
    static TrackFinderFunction makeTrackFinderFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                                       BFieldVariant                                 magneticField);

    /// Type erased track finder function.
    TrackFinderFunction findTracks;

    StatusCode initialize() override;

    StatusCode execute() override;
  };

} // namespace Jug::Reco

#endif
