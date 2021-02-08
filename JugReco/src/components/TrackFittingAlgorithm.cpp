//
#include "TrackFittingAlgorithm.h"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Property.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Definitions/Units.hpp"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugReco/GeometryContainers.hpp"
#include "JugReco/SourceLinks.h"
#include "JugReco/Track.hpp"
#include "JugReco/BField.h"
#include "JugReco/Measurement.hpp"
#include "JugReco/SourceLinks.h"

#include "eicd/TrackerHitCollection.h"

#include <functional>
#include <stdexcept>
#include <vector>
#include <random>
#include <stdexcept>

namespace Jug::Reco {

  using namespace Acts::UnitLiterals;

  TrackFittingAlgorithm::TrackFittingAlgorithm(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc)
  {
    declareProperty("inputSourceLinks", m_inputSourceLinks, "");
    declareProperty("initialTrackParameters", m_initialTrackParameters, "");
    declareProperty("foundTracks", m_foundTracks, "");
    declareProperty("outputTrajectories", m_outputTrajectories, "");
  }

  StatusCode TrackFittingAlgorithm::initialize()
  {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    m_BField   = std::make_shared<Acts::ConstantBField>(Acts::Vector3{0.0, 0.0, m_geoSvc->centralMagneticField()});
    m_fieldctx = BFieldVariant(m_BField);
    // chi2 and #sourclinks per surface cutoffs
    //m_sourcelinkSelectorCfg = {
    //    {Acts::GeometryIdentifier(), {15, 10}},
    //};
    m_trackFittingFunc = makeTrackFittingFunction(m_geoSvc->trackingGeometry(), m_BField);
    return StatusCode::SUCCESS;
  }

  StatusCode TrackFittingAlgorithm::execute()
  {
    // Read input data
    const SourceLinkContainer*      src_links       = m_inputSourceLinks.get();
    const TrackParametersContainer* init_trk_params = m_initialTrackParameters.get();

    // TrajectoryContainer trajectories;
    auto trajectories = m_outputTrajectories.createAndPut();
    trajectories->reserve(init_trk_params->size());

    //// Construct a perigee surface as the target surface
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0., 0., 0.});
    ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("TrackFittingAlgorithm Logger", Acts::Logging::INFO));

    // Perform the track finding for each starting parameter
    // @TODO: use seeds from track seeding algorithm as starting parameter
    for (std::size_t iseed = 0; iseed < init_trk_params->size(); ++iseed) {

      const auto& initialParams = (*init_trk_params)[iseed];
      Acts::PropagatorPlainOptions pOptions;
      pOptions.maxSteps = 10000;

      Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
        m_geoctx, m_fieldctx, m_calibctx,
        Acts::VoidOutlierFinder(), Acts::LoggerWrapper{logger()}, Acts::PropagatorPlainOptions(), &(*pSurface));

      debug() << "Invoke track fitting ...  " << iseed << endmsg;
      auto result = m_trackFittingFunc(*src_links, initialParams, kfOptions);
      debug() << "fitting done." << endmsg;
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
        // trajectories->push_back(SimMultiTrajectory());
      }
    }

    // ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
    return StatusCode::SUCCESS;


    ///////////////////////////
    // acts example

  // Set the KalmanFitter options

  // Perform the fit for each input track
  //std::vector<IndexSourceLink> trackSourceLinks;
  //for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
  //  // The list of hits and the initial start parameters
  //  const auto& protoTrack = protoTracks[itrack];
  //  const auto& initialParams = initialParameters[itrack];

  //  // We can have empty tracks which must give empty fit results so the number
  //  // of entries in input and output containers matches.
  //  if (protoTrack.empty()) {
  //    trajectories.push_back(Trajectories());
  //    ACTS_WARNING("Empty track " << itrack << " found.");
  //    continue;
  //  }

  //  // Clear & reserve the right size
  //  trackSourceLinks.clear();
  //  trackSourceLinks.reserve(protoTrack.size());

  //  // Fill the source links via their indices from the container
  //  for (auto hitIndex : protoTrack) {
  //    auto sourceLink = sourceLinks.nth(hitIndex);
  //    if (sourceLink == sourceLinks.end()) {
  //      ACTS_FATAL("Proto track " << itrack << " contains invalid hit index"
  //                                << hitIndex);
  //      return ProcessCode::ABORT;
  //    }
  //    trackSourceLinks.push_back(*sourceLink);
  //  }

  //  ACTS_DEBUG("Invoke fitter");
  //  auto result = m_cfg.fit(trackSourceLinks, initialParams, kfOptions);
  //  if (result.ok()) {
  //    // Get the fit output object
  //    const auto& fitOutput = result.value();
  //    // The track entry indices container. One element here.
  //    std::vector<size_t> trackTips;
  //    trackTips.reserve(1);
  //    trackTips.emplace_back(fitOutput.trackTip);
  //    // The fitted parameters container. One element (at most) here.
  //    Trajectories::IndexedParameters indexedParams;
  //    if (fitOutput.fittedParameters) {
  //      const auto& params = fitOutput.fittedParameters.value();
  //      ACTS_VERBOSE("Fitted paramemeters for track " << itrack);
  //      ACTS_VERBOSE("  " << params.parameters().transpose());
  //      // Push the fitted parameters to the container
  //      indexedParams.emplace(fitOutput.trackTip, std::move(params));
  //    } else {
  //      ACTS_DEBUG("No fitted paramemeters for track " << itrack);
  //    }
  //    // store the result
  //    trajectories.emplace_back(std::move(fitOutput.fittedStates),
  //                              std::move(trackTips), std::move(indexedParams));
  //  } else {
  //    ACTS_WARNING("Fit failed for track " << itrack << " with error"
  //                                         << result.error());
  //    // Fit failed. Add an empty result so the output container has
  //    // the same number of entries as the input.
  //    trajectories.push_back(Trajectories());
  //  }
  //}


    return StatusCode::SUCCESS;
  }

  DECLARE_COMPONENT(TrackFittingAlgorithm)
} // namespace Jug::Reco
