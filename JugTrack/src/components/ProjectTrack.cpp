#include <algorithm>

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

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/UniqueID.h"
#include "JugBase/BField/DD4hepBField.h"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

// Event Model related classes
#include "eicd/BasicParticleCollection.h"
#include "eicd/TrackerHitCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Trajectories.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

#include "eicd/VectorPolar.h"
#include "eicd/VectorXYZ.h"

#include <cmath>

using BoundTrackParamPtr = std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;

namespace Jug::Reco {

  /** Extrac the particles form fit trajectories.
   *
   * \ingroup tracking
   */
   class ProjectTrack : public GaudiAlgorithm, AlgorithmIDMixin<int32_t> {
   public:
    //DataHandle<eic::RawTrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoriesContainer>     m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::TrackParametersCollection> m_outputTrackParameters{"outputTrackParameters", Gaudi::DataHandle::Writer, this};

    SmartIF<IGeoSvc> m_geoSvc;
    std::shared_ptr<const Acts::TrackingGeometry> m_tGeometry;
    Acts::GeometryContext m_geoContext;
    Acts::MagneticFieldContext m_fieldContext;

   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ProjectTrack(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
        , AlgorithmIDMixin(name, info()) {
          declareProperty("inputTrajectories", m_inputTrajectories,"");
          declareProperty("outputTrackParameters", m_outputTrackParameters, "ACTS Track Parameters");
        }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;

      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }

      m_tGeometry = m_geoSvc->trackingGeometry();

      return StatusCode::SUCCESS;
    }

    BoundTrackParamPtrResult propagateTrack(
        const Acts::BoundTrackParameters& params, const SurfacePtr& targetSurf) {
      
      std::cout << "Propagating final track fit with momentum: " 
            << params.momentum() << " and position " 
            << params.position(m_geoContext)
            << std::endl
            << "track fit phi/eta "
            << atan2(params.momentum()(1), 
               params.momentum()(0)) 
            << " and " 
            << atanh(params.momentum()(2) 
               / params.momentum().norm())
            << std::endl;

      auto magField = std::dynamic_pointer_cast<const Jug::BField::DD4hepBField>(m_geoSvc->getFieldProvider());
      auto bfieldvariant = Jug::BField::BFieldVariant(magField);

      return std::visit([params, targetSurf, this]
            (auto && inputField) -> BoundTrackParamPtrResult {
          using InputMagneticField = 
      typename std::decay_t<decltype(inputField)>::element_type;
          using MagneticField      = Acts::SharedBField<InputMagneticField>;
          using Stepper            = Acts::EigenStepper<MagneticField>;
          using Propagator         = Acts::Propagator<Stepper>;

          MagneticField field(inputField);
          Stepper stepper(field);
          Propagator propagator(stepper);

          Acts::Logging::Level logLevel = Acts::Logging::FATAL;
          // logLevel = Acts::Logging::DEBUG;

          ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("ProjectTrack Logger", logLevel));
          
          Acts::PropagatorOptions<> options(m_geoContext,
              m_fieldContext,
              Acts::LoggerWrapper{logger()});
         
          auto result = propagator.propagate(params, *targetSurf, 
               options);
       
          if(result.ok()) return std::move((*result).endParameters);
         
          return result.error();
       },
         std::move(bfieldvariant));
    }

    StatusCode execute() override {
      // input collection
      const TrajectoriesContainer* trajectories = m_inputTrajectories.get();
      // create output collections
      auto track_pars = m_outputTrackParameters.createAndPut();

      if (msgLevel(MSG::DEBUG)) {
        debug() << std::size(*trajectories) << " trajectories " << endmsg;
      }

      // Loop over the trajectories
      for (size_t itraj = 0; itraj < trajectories->size(); ++itraj) {
        const auto& traj = (*trajectories)[itraj];

        // Get the entry index for the single trajectory
        // The trajectory entry indices and the multiTrajectory
        const auto& mj        = traj.multiTrajectory();
        const auto& trackTips = traj.tips();
        debug() << "# of elements in trackTips " <<trackTips.size() << endmsg;

        if (trackTips.empty()) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "Empty multiTrajectory." << endmsg;
          }
          continue;
        }

        auto& trackTip = trackTips.front();

        // Collect the trajectory summary info
        auto trajState       = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
        int  m_nMeasurements = trajState.nMeasurements;
        int  m_nStates       = trajState.nStates;
        int  m_nCalibrated   = 0;
          
        // get path length at last silicon layer
        float pathlength_at_lastlayer = -9999.;
        mj.visitBackwards(trackTip, [&](auto&& trackstate) {
          // debug() << trackstate.hasPredicted() << endmsg;
          // debug() << trackstate.predicted() << endmsg;
          auto params = trackstate.predicted(); //<< endmsg;
          auto pathlength = trackstate.pathLength();
          auto geoID = trackstate.referenceSurface().geometryId();
          auto volume = geoID.volume();
          auto layer = geoID.layer();
          if (trackstate.hasCalibrated())
          {
            m_nCalibrated++;
          } 

          debug() << "******************************" << endmsg;
          debug() << "predicted variables: \n" << trackstate.predicted() << endmsg;
          debug() << "geoID = " << geoID << endmsg;
          debug() << "volume = " << volume << ", layer = " << layer << endmsg;
          debug() << "pathlength = " << pathlength << endmsg;
          debug() << "hasCalibrated = " << trackstate.hasCalibrated() << endmsg;

          debug() << "pos 0 = " << params[Acts::eFreePos0] << endmsg;
          debug() << "pos 1 = " << params[Acts::eFreePos1] << endmsg;
          debug() << "pos 2 = " << params[Acts::eFreePos2] << endmsg;
          debug() << "******************************" << endmsg;
        });

        // Get the fitted track parameter
        bool m_hasFittedParams = false;
        if (traj.hasTrackParameters(trackTip)) {
          m_hasFittedParams      = true;
          const auto& boundParam = traj.trackParameters(trackTip);
          const auto& parameter  = boundParam.parameters();
          const auto& covariance = *boundParam.covariance();

          if (msgLevel(MSG::DEBUG)) {  
            debug() << "loc 0 = " << parameter[Acts::eBoundLoc0] << endmsg;
            debug() << "loc 1 = " << parameter[Acts::eBoundLoc1] << endmsg;
            debug() << "phi   = " << parameter[Acts::eBoundPhi] << endmsg;
            debug() << "theta = " << parameter[Acts::eBoundTheta] << endmsg;
            debug() << "q/p   = " << parameter[Acts::eBoundQOverP] << endmsg;
            debug() << "p     = " << 1.0 / parameter[Acts::eBoundQOverP] << endmsg;

            debug() << "err phi = " << sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)) << endmsg;
            debug() << "err th  = " << sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)) << endmsg;
            debug() << "err q/p = " << sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)) << endmsg;

            debug() << " chi2 = " << trajState.chi2Sum << endmsg;


            //=================================================
            //               Track projection
            //
            // Reference sPHENIX code: https://github.com/sPHENIX-Collaboration/coresoftware/blob/335e6da4ccacc8374cada993485fe81d82e74a4f/offline/packages/trackreco/PHActsTrackProjection.h
            //=================================================
            // make a reference cylinder to mimic DIRC
            const auto dircRadius = 900;
            const auto eta = 1.1;
            const auto theta = 2. * atan(exp(-eta));
            const auto halfZ = dircRadius / tan(theta);
            
            auto transform = Acts::Transform3::Identity();
      
            std::shared_ptr<Acts::CylinderSurface> dircSurf = 
              Acts::Surface::makeShared<Acts::CylinderSurface>(transform,
                                                               dircRadius,
                                                               halfZ);
            
            // project track parameters to target surface
            auto result = propagateTrack(boundParam, dircSurf);
            if(result.ok())
            {
              auto trackStateParams = std::move(**result);
            }
          }
        }
      }
      
      return StatusCode::SUCCESS;
    }

  };
  DECLARE_COMPONENT(ProjectTrack)

} // namespace Jug::Reco
