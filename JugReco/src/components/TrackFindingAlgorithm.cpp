//
#include "TrackFindingAlgorithm.h"

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

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFinding/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugReco/GeometryContainers.hpp"
#include "JugReco/SourceLinks.h"
#include "JugReco/Track.hpp"
#include "JugReco/BField.h"
#include "JugReco/SourceLinks.h"

#include "eicd/TrackerHitCollection.h"

#include <functional>
#include <stdexcept>
#include <vector>
#include <random>
#include <stdexcept>



namespace Jug::Reco {

  using namespace Acts::UnitLiterals;

  TrackFindingAlgorithm::TrackFindingAlgorithm(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputSourceLinks", m_inputSourceLinks, "");
    declareProperty("inputInitialTrackParameters", m_inputInitialTrackParameters, "");
    declareProperty("outputTrajectories", m_outputTrajectories, "");
  }

  StatusCode TrackFindingAlgorithm::initialize() {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    m_BField                = std::make_shared<Acts::ConstantBField>(Acts::Vector3D{0.0, 0.0, m_geoSvc->centralMagneticField()});
    m_fieldctx              = BFieldVariant(m_BField);

    // chi2 and #sourclinks per surface cutoffs
    m_sourcelinkSelectorCfg = { {Acts::GeometryIdentifier(), {15, 10}},
    };

    findTracks = TrackFindingAlgorithm::makeTrackFinderFunction(m_geoSvc->trackingGeometry(), m_BField);

    return StatusCode::SUCCESS;
  }

  StatusCode TrackFindingAlgorithm::execute() {
    // Read input data
    const SourceLinkContainer*      src_links       = m_inputSourceLinks.get();
    const TrackParametersContainer* init_trk_params = m_inputInitialTrackParameters.get();
    // const auto sourceLinks       = ctx.eventStore.get<SourceLinkContainer>(m_cfg.inputSourceLinks);
    // const auto initialParameters = ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputInitialTrackParameters);

    //// Prepare the output data with MultiTrajectory
    // TrajectoryContainer trajectories;
    auto trajectories = m_outputTrajectories.createAndPut();
    trajectories->reserve(init_trk_params->size());

    //// Construct a perigee surface as the target surface
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3D{0., 0., 0.});

    ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("TrackFindingAlgorithm Logger", Acts::Logging::INFO));

    // Perform the track finding for each starting parameter
    // @TODO: use seeds from track seeding algorithm as starting parameter
    for (std::size_t iseed = 0; iseed < init_trk_params->size(); ++iseed) {
      const auto& initialParams = (*init_trk_params)[iseed];

      Acts::PropagatorPlainOptions pOptions;
      pOptions.maxSteps = 10000;
      // Set the CombinatorialKalmanFilter options
      CKFOptions ckfOptions(m_geoctx, m_fieldctx, m_calibctx, m_sourcelinkSelectorCfg,
                                      Acts::LoggerWrapper{logger()}, pOptions,
                                      &(*pSurface));
      // TrackFindingAlgorithm::CKFOptions ckfOptions(ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
      //                                                 m_cfg.sourcelinkSelectorCfg, Acts::LoggerWrapper{logger()},
      //                                                 &(*pSurface));

      debug() << "Invoke track finding seeded by truth particle " << iseed << endmsg;
      auto result = findTracks(*src_links, initialParams, ckfOptions);
      debug() << "finding done." << endmsg;
      if (result.ok()) {
        // Get the track finding output object
        const auto& trackFindingOutput = result.value();
        // Create a SimMultiTrajectory
        trajectories->emplace_back(std::move(trackFindingOutput.fittedStates), std::move(trackFindingOutput.trackTips),
                                   std::move(trackFindingOutput.fittedParameters));
      } else {
        debug() << "Track finding failed for truth seed " << iseed << endmsg;
        ACTS_WARNING("Track finding failed for truth seed " << iseed << " with error" << result.error());
        // Track finding failed, but still create an empty SimMultiTrajectory
        //trajectories->push_back(SimMultiTrajectory());
      }
    }

    // ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
    return StatusCode::SUCCESS;
  }

  DECLARE_COMPONENT(TrackFindingAlgorithm)
} // namespace Jug::Reco

namespace {
  template <typename TrackFinder>
  struct TrackFinderFunctionImpl {
    TrackFinder trackFinder;

    TrackFinderFunctionImpl(TrackFinder&& f) : trackFinder(std::move(f)) {}

    Jug::Reco::TrackFindingAlgorithm::TrackFinderResult
    operator()(const Jug::SourceLinkContainer& sourceLinks, const Jug::TrackParameters& initialParameters,
               const Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>& options) const {
      return trackFinder.findTracks(sourceLinks, initialParameters, options);
    };
  };
} // namespace

namespace Jug::Reco {

  TrackFindingAlgorithm::TrackFinderFunction
  TrackFindingAlgorithm::makeTrackFinderFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                                 BFieldVariant                                 magneticField)
  {
    using Updater  = Acts::GainMatrixUpdater;
    using Smoother = Acts::GainMatrixSmoother;

    // unpack the magnetic field variant and instantiate the corresponding track
    // finder.
    return std::visit(
        [trackingGeometry](auto&& inputField) -> TrackFinderFunction {
          // each entry in the variant is already a shared_ptr
          // need ::element_type to get the real magnetic field type
          using InputMagneticField = typename std::decay_t<decltype(inputField)>::element_type;
          using MagneticField      = Acts::SharedBField<InputMagneticField>;
          using Stepper            = Acts::EigenStepper<MagneticField>;
          using Navigator          = Acts::Navigator;
          using Propagator         = Acts::Propagator<Stepper, Navigator>;
          using SourceLinkSelector = Acts::CKFSourceLinkSelector;
          using CKF                = Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother, SourceLinkSelector>;

          std::cout << " finding ...\n";
          // construct all components for the track finder
          MagneticField field(std::move(inputField));
          Stepper       stepper(std::move(field));
          Navigator     navigator(trackingGeometry);
          navigator.resolvePassive   = false;
          navigator.resolveMaterial  = true;
          navigator.resolveSensitive = true;
          std::cout << " propagator\n" ;
          Propagator propagator(std::move(stepper), std::move(navigator));
          CKF        trackFinder(std::move(propagator));

          // build the track finder functions. owns the track finder object.
          return TrackFinderFunctionImpl<CKF>(std::move(trackFinder));
        },
        std::move(magneticField));
  }
} // namespace Jug::Reco

