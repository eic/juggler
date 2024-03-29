// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

#include "CKFTracking.h"

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

#include "Acts/TrackFinding/SourceLinkAccessorConcept.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Definitions/Units.hpp"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>
#include "JugTrack/IActsGeoSvc.h"
#include "JugTrack/DD4hepBField.h"

#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "edm4eic/TrackerHitCollection.h"

#include <functional>
#include <stdexcept>
#include <system_error>
#include <type_traits>
#include <vector>
#include <random>
#include <stdexcept>

template<> struct fmt::formatter<std::error_code> : fmt::ostream_formatter {};

static const std::map<int, Acts::Logging::Level> s_msgMap = {
    {MSG::DEBUG, Acts::Logging::DEBUG},
    {MSG::VERBOSE, Acts::Logging::VERBOSE},
    {MSG::INFO, Acts::Logging::INFO},
    {MSG::WARNING, Acts::Logging::WARNING},
    {MSG::FATAL, Acts::Logging::FATAL},
    {MSG::ERROR, Acts::Logging::ERROR},
};

namespace Jug::Reco {

  using namespace Acts::UnitLiterals;

  CKFTracking::CKFTracking(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc)
  {
    declareProperty("inputSourceLinks", m_inputSourceLinks, "");
    declareProperty("inputMeasurements", m_inputMeasurements, "");
    declareProperty("inputInitialTrackParameters", m_inputInitialTrackParameters, "");
    declareProperty("outputTracks", m_outputTracks, "");
    declareProperty("outputTrajectories", m_outputTrajectories, "");
  }

  StatusCode CKFTracking::initialize()
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
    m_actsGeoSvc = service("ActsGeoSvc");
    if (!m_actsGeoSvc) {
      error() << "Unable to locate ACTS Geometry Service. "
              << "Make sure you have ActsGeoSvc in the right place in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }

    m_BField   = std::dynamic_pointer_cast<const Jug::BField::DD4hepBField>(m_actsGeoSvc->getFieldProvider());
    m_fieldctx = Jug::BField::BFieldVariant(m_BField);

    // eta bins, chi2 and #sourclinks per surface cutoffs
    m_sourcelinkSelectorCfg = {
        {Acts::GeometryIdentifier(),
            {m_etaBins, m_chi2CutOff,
                {m_numMeasurementsCutOff.begin(), m_numMeasurementsCutOff.end()}
            }
        },
    };
    m_trackFinderFunc = CKFTracking::makeCKFTrackingFunction(m_actsGeoSvc->trackingGeometry(), m_BField);
    auto im = s_msgMap.find(msgLevel());
    if (im != s_msgMap.end()) {
        m_actsLoggingLevel = im->second;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode CKFTracking::execute()
  {
    // Read input data
    const auto* const src_links       = m_inputSourceLinks.get();
    const auto* const init_trk_params = m_inputInitialTrackParameters.get();
    const auto* const measurements    = m_inputMeasurements.get();

    //// Prepare the output data with MultiTrajectory
    auto* acts_trajectories = m_outputTrajectories.createAndPut();
    acts_trajectories->reserve(init_trk_params->size());

    //// Construct a perigee surface as the target surface
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0., 0., 0.});

    ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("CKFTracking Logger", m_actsLoggingLevel));

    Acts::PropagatorPlainOptions pOptions;
    pOptions.maxSteps = 10000;

    ActsExamples::PassThroughCalibrator pcalibrator;
    ActsExamples::MeasurementCalibratorAdapter calibrator(pcalibrator, *measurements);
    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;
    Acts::MeasurementSelector measSel{m_sourcelinkSelectorCfg};

    Acts::CombinatorialKalmanFilterExtensions<Acts::VectorMultiTrajectory>
        extensions;
    extensions.calibrator.connect<
        &ActsExamples::MeasurementCalibratorAdapter::calibrate>(
        &calibrator);
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &kfUpdater);
    extensions.smoother.connect<
        &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
        &kfSmoother);
    extensions.measurementSelector
        .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
            &measSel);

    ActsExamples::IndexSourceLinkAccessor slAccessor;
    slAccessor.container = src_links;
    Acts::SourceLinkAccessorDelegate<ActsExamples::IndexSourceLinkAccessor::Iterator>
        slAccessorDelegate;
    slAccessorDelegate.connect<&ActsExamples::IndexSourceLinkAccessor::range>(&slAccessor);

    // Set the CombinatorialKalmanFilter options
    CKFTracking::TrackFinderOptions options(
        m_geoctx, m_fieldctx, m_calibctx, slAccessorDelegate,
        extensions, pOptions, &(*pSurface));

    // Create track container
    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    ActsExamples::TrackContainer tracks(trackContainer, trackStateContainer);

    // Add seed number column
    tracks.addColumn<unsigned int>("seed");
    Acts::TrackAccessor<unsigned int> seedNumber("seed");

    // Loop over seeds
    for (std::size_t iseed = 0; iseed < init_trk_params->size(); ++iseed) {
        auto result =
            (*m_trackFinderFunc)(init_trk_params->at(iseed), options, tracks);

        if (!result.ok()) {
            debug() << fmt::format("Track finding failed for seed {} with error {}", iseed, result.error()) << endmsg;
            continue;
        }

        // Set seed number for all found tracks
        auto& tracksForSeed = result.value();
        for (auto& track : tracksForSeed) {
            seedNumber(track) = iseed;
        }
    }

    // Move track states and track container to const containers
    // NOTE Using the non-const containers leads to references to
    // implicitly converted temporaries inside the Trajectories.
    auto constTrackStateContainer =
        std::make_shared<Acts::ConstVectorMultiTrajectory>(
            std::move(*trackStateContainer));

    auto constTrackContainer =
        std::make_shared<Acts::ConstVectorTrackContainer>(
            std::move(*trackContainer));

    auto* constTracks = m_outputTracks.put(
        std::make_unique<ActsExamples::ConstTrackContainer>(constTrackContainer, constTrackStateContainer)
    );

    // Seed number column accessor
    const Acts::ConstTrackAccessor<unsigned int> constSeedNumber("seed");


    // Prepare the output data with MultiTrajectory, per seed
    ActsExamples::Trajectories::IndexedParameters parameters;
    std::vector<Acts::MultiTrajectoryTraits::IndexType> tips;

    std::optional<unsigned int> lastSeed;
    for (const auto& track : *constTracks) {
      if (!lastSeed) {
        lastSeed = constSeedNumber(track);
      }

      if (constSeedNumber(track) != lastSeed.value()) {
        // make copies and clear vectors
        acts_trajectories->emplace_back(
          constTracks->trackStateContainer(),
          tips, parameters);

        tips.clear();
        parameters.clear();
      }

      lastSeed = constSeedNumber(track);

      tips.push_back(track.tipIndex());
      parameters.emplace(
          std::pair{track.tipIndex(),
                    ActsExamples::TrackParameters{track.referenceSurface().getSharedPtr(),
                                                  track.parameters(), track.covariance(),
                                                  track.particleHypothesis()}});
    }

    if (tips.empty()) {
      info() << fmt::format("Last trajectory is empty") << endmsg;
    }

    // last entry: move vectors
    acts_trajectories->emplace_back(
      constTracks->trackStateContainer(),
      std::move(tips), std::move(parameters));

    return StatusCode::SUCCESS;
  }

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(CKFTracking)
} // namespace Jug::Reco
