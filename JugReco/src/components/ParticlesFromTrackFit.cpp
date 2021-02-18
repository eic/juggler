#include <algorithm>

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

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
 
// Event Model related classes
#include "eicd/ParticleCollection.h"
#include "eicd/TrackerHitCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "JugReco/IndexSourceLink.hpp"
#include "JugReco/Track.hpp"

#include "Acts/Utilities/Helpers.hpp"


namespace Jug {
  namespace Reco {
  
    /** Ultra-fast silicon detector digitization.
     *
     */
   class ParticlesFromTrackFit : public GaudiAlgorithm {
   public:
    //DataHandle<eic::RawTrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoryContainer>     m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer, this};
    DataHandle<eic::TrackParametersCollection> m_outputTrackParameters{"outputTrackParameters", Gaudi::DataHandle::Writer, this};

   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ParticlesFromTrackFit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputTrajectories", m_inputTrajectories,"");
          declareProperty("outputParticles", m_outputParticles, "");
          declareProperty("outputTrackParameters", m_outputTrackParameters, "ACTS Track Parameters");
        }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const TrajectoryContainer* trajectories = m_inputTrajectories.get();
      // create output collections
      auto rec_parts = m_outputParticles.createAndPut();
      auto track_pars = m_outputTrackParameters.createAndPut();

      debug() << std::size(*trajectories) << " trajectories " << endmsg;

      for(const auto& traj : *trajectories) {
        //traj.trajectory().first
        const auto& [trackTips, mj] = traj.trajectory();
        if (trackTips.empty()) {
          debug() << "Empty multiTrajectory." << endmsg;
          continue;
        }

        // Get the entry index for the single trajectory
        auto& trackTip = trackTips.front();

        // Collect the trajectory summary info
        auto trajState       = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
        int m_nMeasurements = trajState.nMeasurements;
        int m_nStates = trajState.nStates;

        // Get the fitted track parameter
        bool  m_hasFittedParams = false;
        if (traj.hasTrackParameters(trackTip)) {
          m_hasFittedParams = true;
          const auto& boundParam = traj.trackParameters(trackTip);
          const auto& parameter = boundParam.parameters();
          const auto& covariance = *boundParam.covariance();
          debug() << "loc 0 = " << parameter[Acts::eBoundLoc0]   << endmsg;
          debug() << "loc 1 = " << parameter[Acts::eBoundLoc1]   << endmsg;
          debug() << "phi   = " << parameter[Acts::eBoundPhi]    << endmsg;
          debug() << "theta = " << parameter[Acts::eBoundTheta]  << endmsg;
          debug() << "q/p   = " << parameter[Acts::eBoundQOverP] << endmsg;
          debug() << "p     = " << 1.0/parameter[Acts::eBoundQOverP] << endmsg;

          debug() << "err phi = " << sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi))      << endmsg;
          debug() << "err th  = " << sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta))  << endmsg;
          debug() << "err q/p = " << sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP))<< endmsg;

          debug() << " chi2 = " << trajState.chi2Sum << endmsg;

          eic::TrackParameters pars({
            parameter[Acts::eBoundLoc0], parameter[Acts::eBoundLoc1], parameter[Acts::eBoundPhi],
            parameter[Acts::eBoundTheta], parameter[Acts::eBoundQOverP],parameter[Acts::eBoundTime],
            sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)), sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)),
            sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)), sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
            sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)),sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime))});
          track_pars->push_back(pars);

          //m_ePHI_fit = parameter[Acts::eBoundPhi];
          //m_eTHETA_fit = parameter[Acts::eBoundTheta];
          //m_eQOP_fit = parameter[Acts::eBoundQOverP];
          //m_eT_fit = parameter[Acts::eBoundTime];
          //m_err_eLOC0_fit = sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0));
          //m_err_eLOC1_fit = sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1));
          //m_err_ePHI_fit = sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi));
          //m_err_eTHETA_fit = sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta));
          //m_err_eQOP_fit = sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP));
          //m_err_eT_fit = sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime));
        }

        auto tsize = traj.trajectory().first.size();
        debug() << "# fitted parameters : " << tsize << endmsg;
        if(tsize == 0 ) continue;

        traj.trajectory().second.visitBackwards(tsize-1, [&](auto&& trackstate) {
          //debug() << trackstate.hasPredicted() << endmsg;
          //debug() << trackstate.predicted() << endmsg;
          auto params = trackstate.predicted() ;//<< endmsg;

          double p0 = (1.0 / params[Acts::eBoundQOverP]) / Acts::UnitConstants::GeV;
          debug() << "track predicted p = " << p0 << " GeV" << endmsg;
          if ( std::abs(p0) > 500) {
            debug() << "skipping" << endmsg;
            return;
          }

          eic::Particle p({params[Acts::eBoundPhi], params[Acts::eBoundTheta], 1.0 / std::abs(params[Acts::eBoundQOverP]), 0.000511},
                          {0.0, 0.0, 0.0, params[Acts::eBoundTime]},
                          (long long)11 * params[Acts::eBoundQOverP] / std::abs(params[Acts::eBoundQOverP]), 0);
          //debug() << p << endmsg;
          rec_parts->push_back(p);
        });

      }
      return StatusCode::SUCCESS;
    }

  };
  DECLARE_COMPONENT(ParticlesFromTrackFit)

  } // namespace Examples
} // namespace Gaudi

