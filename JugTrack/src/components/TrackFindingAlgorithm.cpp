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

#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
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
#include "JugBase/BField/DD4hepBField.h"

#include "JugTrack/GeometryContainers.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/Index.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"


#include "eicd/TrackerHitCollection.h"

#include <functional>
#include <stdexcept>
#include <vector>
#include <random>
#include <stdexcept>


static const std::map<int, Acts::Logging::Level> _msgMap = {
    {MSG::DEBUG, Acts::Logging::DEBUG},
    {MSG::VERBOSE, Acts::Logging::VERBOSE},
    {MSG::INFO, Acts::Logging::INFO},
    {MSG::WARNING, Acts::Logging::WARNING},
    {MSG::FATAL, Acts::Logging::FATAL},
    {MSG::ERROR, Acts::Logging::ERROR},
};

namespace Jug::Reco {

  using namespace Acts::UnitLiterals;

  TrackFindingAlgorithm::TrackFindingAlgorithm(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc)
  {
    declareProperty("inputSourceLinks", m_inputSourceLinks, "");
    declareProperty("inputMeasurements", m_inputMeasurements, "");
    declareProperty("inputInitialTrackParameters", m_inputInitialTrackParameters, "");
    declareProperty("outputTrajectories", m_outputTrajectories, "");
  }

  StatusCode TrackFindingAlgorithm::initialize()
  {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    //m_BField   = m_geoSvc->getFieldProvider();//std::make_shared<Acts::ConstantBField>(Acts::Vector3{0.0, 0.0, m_geoSvc->centralMagneticField()});
    m_BField   = std::dynamic_pointer_cast<const Jug::BField::DD4hepBField>(m_geoSvc->getFieldProvider());
    m_fieldctx = Jug::BField::BFieldVariant(m_BField);


    for(int z : {0,1000,2000,4000,5899}){
      auto b =  m_BField->getField({0.0,0.0,double(z)})/(Acts::UnitConstants::T);
      debug() << "B(z=" << z << " mm) = " << b.transpose()  << " T"   << endmsg;
    }



    // chi2 and #sourclinks per surface cutoffs
    m_sourcelinkSelectorCfg = {
        {Acts::GeometryIdentifier(), {15, 10}},
    };
    m_trackFinderFunc = TrackFindingAlgorithm::makeTrackFinderFunction(m_geoSvc->trackingGeometry(), m_BField);
    auto im = _msgMap.find(msgLevel());
    if (im != _msgMap.end()) {
        m_actsLoggingLevel = im->second;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode TrackFindingAlgorithm::execute()
  {
    // Read input data
    const IndexSourceLinkContainer* src_links       = m_inputSourceLinks.get();
    const TrackParametersContainer* init_trk_params = m_inputInitialTrackParameters.get();
    const MeasurementContainer*     measurements    = m_inputMeasurements.get();

    //// Prepare the output data with MultiTrajectory
    // TrajectoryContainer trajectories;
    auto trajectories = m_outputTrajectories.createAndPut();
    trajectories->reserve(init_trk_params->size());

    //// Construct a perigee surface as the target surface
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0., 0., 0.});

    ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("TrackFindingAlgorithm Logger", m_actsLoggingLevel));

    // Perform the track finding for each starting parameter
    // @TODO: use seeds from track seeding algorithm as starting parameter
    // for (std::size_t iseed = 0; iseed < init_trk_params->size(); ++iseed) {
    //  const auto& initialParams = (*init_trk_params)[iseed];

    Acts::PropagatorPlainOptions pOptions;
    pOptions.maxSteps = 10000;

    // Set the CombinatorialKalmanFilter options
    TrackFindingAlgorithm::TrackFinderOptions options(
        m_geoctx, m_fieldctx, m_calibctx, MeasurementCalibrator(*measurements),
        Acts::MeasurementSelector(m_sourcelinkSelectorCfg), Acts::LoggerWrapper{logger()}, pOptions, &(*pSurface));

    auto results = m_trackFinderFunc(*src_links, *init_trk_params, options);

    for (std::size_t iseed = 0; iseed < init_trk_params->size(); ++iseed) {

      auto& result = results[iseed];

      if (result.ok()) {
        // Get the track finding output object
        const auto& trackFindingOutput = result.value();
        // Create a SimMultiTrajectory
        trajectories->emplace_back(std::move(trackFindingOutput.fittedStates), 
                                   std::move(trackFindingOutput.trackTips),
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
  }

  DECLARE_COMPONENT(TrackFindingAlgorithm)
} // namespace Jug::Reco

