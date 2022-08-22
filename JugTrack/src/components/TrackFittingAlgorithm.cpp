// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

//
#include "TrackFittingAlgorithm.h"

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
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Definitions/Units.hpp"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/BField/DD4hepBField.h"

#include "JugTrack/GeometryContainers.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Measurement.hpp"

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
    declareProperty("inputMeasurements", m_inputMeasurements, "");
    declareProperty("inputProtoTracks", m_inputProtoTracks, "");
    declareProperty("foundTracks", m_foundTracks, "");
    declareProperty("outputTrajectories", m_outputTrajectories, "");
  }

  StatusCode TrackFittingAlgorithm::initialize()
  {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    m_BField   = std::dynamic_pointer_cast<const Jug::BField::DD4hepBField>(m_geoSvc->getFieldProvider());
    m_fieldctx = Jug::BField::BFieldVariant(m_BField);

    m_trackFittingFunc = makeTrackFittingFunction(m_geoSvc->trackingGeometry(), m_BField);
    return StatusCode::SUCCESS;
  }

  StatusCode TrackFittingAlgorithm::execute()
  {
    // Read input data
    const auto* const sourceLinks       = m_inputSourceLinks.get();
    const auto* const initialParameters = m_initialTrackParameters.get();
    const auto* const measurements      = m_inputMeasurements.get();
    const auto* const protoTracks       = m_inputProtoTracks.get();
    ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("TrackFittingAlgorithm Logger", Acts::Logging::INFO));

    // Consistency cross checks
    if (protoTracks->size() != initialParameters->size()) {
      ACTS_FATAL("Inconsistent number of proto tracks and initial parameters");
      return StatusCode::FAILURE;
    }

    // TrajectoryContainer trajectories;
    auto* trajectories = m_outputTrajectories.createAndPut();
    trajectories->reserve(initialParameters->size());

    // Construct a perigee surface as the target surface
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0., 0., 0.});

    Acts::PropagatorPlainOptions pOptions;
    pOptions.maxSteps = 10000;
    // kfOptions.multipleScattering = m_cfg.multipleScattering;
    // kfOptions.energyLoss         = m_cfg.energyLoss;

    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
    MeasurementCalibrator calibrator{*measurements};
    extensions.calibrator.connect<&MeasurementCalibrator::calibrate>(&calibrator);
    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;
    extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &kfUpdater);
    extensions.smoother.connect<&Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
        &kfSmoother);

    Acts::KalmanFitterOptions kfOptions(
        m_geoctx, m_fieldctx, m_calibctx, extensions,
        Acts::LoggerWrapper{logger()}, Acts::PropagatorPlainOptions(),
        &(*pSurface));

    // used for processing the data
    std::vector<IndexSourceLink>      trackSourceLinks;
    std::vector<const Acts::Surface*> surfSequence;

    if (msgLevel(MSG::DEBUG)) {
      debug() << "initialParams size:  " << initialParameters->size() << endmsg;
      debug() << "measurements size:  " << measurements->size() << endmsg;
    }

    // Perform the track finding for each starting parameter
    // @TODO: use seeds from track seeding algorithm as starting parameter
    // initial track params and proto tracks might likely have the same size.
    //for (std::size_t iseed = 0; iseed < init_trk_params->size(); ++iseed) {
    for (std::size_t itrack = 0; itrack < (*protoTracks).size(); ++itrack) {

      const auto& protoTrack    = (*protoTracks)[itrack];
      const auto& initialParams = (*initialParameters)[itrack];

      if (msgLevel(MSG::DEBUG)) {
        debug() << "protoTrack size:  " << protoTrack.size() << endmsg;
        debug() << "sourceLinks size:  " << sourceLinks->size() << endmsg;
      }

      trackSourceLinks.clear();
      trackSourceLinks.reserve(protoTrack.size());

      for (auto hitIndex : protoTrack) {
        if (auto it = sourceLinks->nth(hitIndex); it != sourceLinks->end()) {
          const IndexSourceLink& sourceLink = *it;
          trackSourceLinks.push_back(std::cref(sourceLink));
        } else {
          ACTS_FATAL("Proto track " << itrack << " contains invalid hit index"
                                    << hitIndex);
          return StatusCode::FAILURE;
        }
      }

      if (msgLevel(MSG::DEBUG)) {
        debug() << "Invoke track fitting ...  " << itrack << endmsg;
      }
      auto result = fitTrack(trackSourceLinks, initialParams, kfOptions);
      if (msgLevel(MSG::DEBUG)) {
        debug() << "fitting done." << endmsg;
      }
      if (result.ok())
      {
        // Get the fit output object
        auto& fitOutput = result.value();
        // The track entry indices container. One element here.
        std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
        trackTips.reserve(1);
        trackTips.emplace_back(fitOutput.lastMeasurementIndex);
        // The fitted parameters container. One element (at most) here.
        Trajectories::IndexedParameters indexedParams;
        if (fitOutput.fittedParameters) {
          const auto& params = fitOutput.fittedParameters.value();
          ACTS_VERBOSE("Fitted paramemeters for track " << itrack);
          ACTS_VERBOSE("  " << params.parameters().transpose());
          // Push the fitted parameters to the container
          indexedParams.emplace(fitOutput.lastMeasurementIndex, params);
        } else {
          ACTS_DEBUG("No fitted paramemeters for track " << itrack);
        }
        // store the result
        trajectories->emplace_back(std::move(fitOutput.fittedStates), std::move(trackTips),
                                  std::move(indexedParams));
      } else {
        ACTS_WARNING("Fit failed for track " << itrack << " with error" << result.error());
        // Fit failed. Add an empty result so the output container has
        // the same number of entries as the input.
        trajectories->push_back(Trajectories());
      }
    }

    return StatusCode::SUCCESS;
  }

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(TrackFittingAlgorithm)

} // namespace Jug::Reco
