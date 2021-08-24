#ifndef JUGGLER_JUGRECO_TrackFittingAlgorithm_HH
#define JUGGLER_JUGRECO_TrackFittingAlgorithm_HH 1

#include <functional>
#include <stdexcept>
#include <vector>
#include <random>
#include <stdexcept>

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/BField/DD4hepBField.h"
#include "JugTrack/GeometryContainers.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/Trajectories.hpp"
#include "JugTrack/ProtoTrack.hpp"

#include "eicd/TrackerHitCollection.h"

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
  class TrackFittingAlgorithm : public GaudiAlgorithm {
  public:
    /// Track fitter function that takes input measurements, initial trackstate
    /// and fitter options and returns some track-fitter-specific result.
    using TrackFitterOptions =
        Acts::KalmanFitterOptions<MeasurementCalibrator, Acts::VoidOutlierFinder>;

    //using TrackFinderResult = Acts::Result<Acts::CombinatorialKalmanFilterResult<SourceLink>>;
    using FitterResult = Acts::Result<Acts::KalmanFitterResult<IndexSourceLink>>;

    /// Fit function that takes input measurements, initial trackstate and fitter
    using FitterFunction = std::function<FitterResult(
      const std::vector<IndexSourceLink>&, const TrackParameters&, const TrackFitterOptions&)>;

  public:
    DataHandle<IndexSourceLinkContainer> m_inputSourceLinks{"inputSourceLinks", Gaudi::DataHandle::Reader, this};
    DataHandle<TrackParametersContainer> m_initialTrackParameters{"initialTrackParameters", Gaudi::DataHandle::Reader, this};
    DataHandle<MeasurementContainer>     m_inputMeasurements{"inputMeasurements", Gaudi::DataHandle::Reader, this};
    DataHandle<ProtoTrackContainer>      m_inputProtoTracks{"inputProtoTracks", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoriesContainer>    m_foundTracks{"foundTracks", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoriesContainer>    m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};

    FitterFunction                        m_trackFittingFunc;
    SmartIF<IGeoSvc>                      m_geoSvc;
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
    // BFieldVariant                                 magneticField);

    StatusCode initialize() override;

    StatusCode execute() override;
   private:
    /// Helper function to call correct FitterFunction
    FitterResult fitTrack(
        const std::vector<IndexSourceLink>& sourceLinks,
        const TrackParameters& initialParameters,
        const TrackFitterOptions& options
        ) const;
        //, const std::vector<const Acts::Surface*>& surfSequence) const;
  };

  inline TrackFittingAlgorithm::FitterResult
  TrackFittingAlgorithm::fitTrack(const std::vector<IndexSourceLink>& sourceLinks,
                                  const TrackParameters&              initialParameters,
                                  const TrackFitterOptions&           options) const
  {
    // const std::vector<const Acts::Surface*>& surfSequence) const
    // if (m_cfg.directNavigation) {
    //  return m_cfg.dFit(sourceLinks, initialParameters, options, surfSequence);
    //}
    return m_trackFittingFunc(sourceLinks, initialParameters, options);
  }

} // namespace Jug::Reco

#endif
