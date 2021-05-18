#ifndef JUGGLER_JUGRECO_TrackFittingAlgorithm_HH
#define JUGGLER_JUGRECO_TrackFittingAlgorithm_HH

#include "JugTrack/GeometryContainers.hpp"

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
//#include "Acts/Definitions/Common.hpp"
//#include "Acts/Utilities/Helpers.hpp"
//#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <stdexcept>
#include <vector>

#include "JugTrack/SourceLinks.h"
#include "JugTrack/Track.hpp"
#include "JugTrack/BField.h"
#include "JugTrack/Measurement.hpp"

#include "eicd/TrackerHitCollection.h"

//#include "Acts/Surfaces/PerigeeSurface.hpp"

#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Definitions/Common.hpp"

//#include "Acts/Fitter/GainMatrixSmoother.hpp"
//#include "Acts/Fitter/GainMatrixUpdater.hpp"
//#include "Acts/MagneticField/ConstantBField.hpp"
//#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
//#include "Acts/MagneticField/SharedBField.hpp"
//#include "Acts/Propagator/EigenStepper.hpp"
//#include "Acts/Propagator/Navigator.hpp"
//#include "Acts/Propagator/Propagator.hpp"
//#include "Acts/Definitions/Units.hpp"

#include <random>
#include <stdexcept>

namespace Jug::Reco {

  class TrackFittingAlgorithm : public GaudiAlgorithm {
  public:
    //using TrackFinderResult = Acts::Result<Acts::CombinatorialKalmanFilterResult<SourceLink>>;
    using FitterResult = Acts::Result<Acts::KalmanFitterResult<SourceLink>>;
    /// Fit function that takes input measurements, initial trackstate and fitter

    using FitterFunction = std::function<FitterResult(
      const std::vector<SourceLink>&, const TrackParameters&,
      const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>&)>;


  /// Track fitter function that takes input measurements, initial trackstate
  /// and fitter options and returns some track-fitter-specific result.
  using TrackFitterOptions =
      Acts::KalmanFitterOptions< Acts::VoidOutlierFinder>;
  using TrackFitterResult =
      Acts::Result<Acts::KalmanFitterResult<SourceLink>>;
  using TrackFitterFunction = std::function<TrackFitterResult(
      const std::vector<SourceLink>&, const TrackParameters&,
      const TrackFitterOptions&)>;

  public:
    DataHandle<SourceLinkContainer>      m_inputSourceLinks{"inputSourceLinks", Gaudi::DataHandle::Reader, this};
    DataHandle<TrackParametersContainer> m_initialTrackParameters{"initialTrackParameters", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoryContainer>      m_foundTracks{"foundTracks", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoryContainer>      m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};

    FitterFunction                  m_trackFittingFunc;
    SmartIF<IGeoSvc>                      m_geoSvc;
    std::shared_ptr<Acts::ConstantBField> m_BField = nullptr;
    Acts::GeometryContext                 m_geoctx;
    Acts::CalibrationContext              m_calibctx;
    Acts::MagneticFieldContext            m_fieldctx;

    //Acts::CKFSourceLinkSelector::Config m_sourcelinkSelectorCfg;

    TrackFittingAlgorithm(const std::string& name, ISvcLocator* svcLoc);

    /** Create the track finder function implementation.
     *  The magnetic field is intentionally given by-value since the variant
     *  contains shared_ptr anyways.
     */
    static FitterFunction makeTrackFittingFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                                        BFieldVariant                                 magneticField);


    StatusCode initialize() override;

    StatusCode execute() override;
  };

} // namespace Jug::Reco

#endif
