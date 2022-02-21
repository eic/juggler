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

#include "eicd/vector_utils.h"

#include <cmath>

namespace Jug::Reco {

  /** Extrac the particles form fit trajectories.
   *
   * \ingroup tracking
   */
   class ParticlesFromTrackFit : public GaudiAlgorithm {
   public:
    //DataHandle<eic::RawTrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoriesContainer>     m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::BasicParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer, this};
    DataHandle<eic::TrackParametersCollection> m_outputTrackParameters{"outputTrackParameters", Gaudi::DataHandle::Writer, this};

   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ParticlesFromTrackFit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputTrajectories", m_inputTrajectories,"");
          declareProperty("outputParticles", m_outputParticles, "");
          declareProperty("outputTrackParameters", m_outputTrackParameters, "Acts Track Parameters");
        }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const TrajectoriesContainer* trajectories = m_inputTrajectories.get();
      // create output collections
      auto rec_parts = m_outputParticles.createAndPut();
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
          if (trackTips.empty()) {
            if (msgLevel(MSG::DEBUG)) {
              debug() << "Empty multiTrajectory." << endmsg;
            }
            continue;
          }

          auto& trackTip = trackTips.front();

          // Collect the trajectory summary info
          auto trajState       = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
          //int  m_nMeasurements = trajState.nMeasurements;
          //int  m_nStates       = trajState.nStates;

          // Get the fitted track parameter
          //
          bool hasFittedParams = false;
          if (traj.hasTrackParameters(trackTip)) {
            hasFittedParams      = true;
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
            }

            const eic::CovXYZ covMomentum {
              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
              static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundPhi)),
              static_cast<float>(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)),
              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundPhi)),
              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundQOverP)),
              static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundQOverP))};
            const eic::CovXY covPos {
              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)),
              static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)),
              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc1))};
            const float timeError{sqrt(static_cast<float>(covariance(Acts::eBoundTime, Acts::eBoundTime)))};

                eic::TrackParameters pars{0, // type: track head --> 0
                                          {parameter[Acts::eBoundLoc0], parameter[Acts::eBoundLoc1]},
                                          covPos,
                                          static_cast<float>(parameter[Acts::eBoundTheta]),
                                          static_cast<float>(parameter[Acts::eBoundPhi]),
                                          static_cast<float>(parameter[Acts::eBoundQOverP]),
                                          covMomentum,
                                          static_cast<float>(parameter[Acts::eBoundTime]),
                                          timeError,
                                          static_cast<float>(boundParam.charge())};
            track_pars->push_back(pars);
          }

          auto tsize = trackTips.size();
          if (msgLevel(MSG::DEBUG)) {
            debug() << "# fitted parameters : " << tsize << endmsg;
          }
          if (tsize == 0)
            continue;

          mj.visitBackwards(tsize - 1, [&](auto&& trackstate) {
            // debug() << trackstate.hasPredicted() << endmsg;
            // debug() << trackstate.predicted() << endmsg;
            auto params = trackstate.predicted(); //<< endmsg;

            double p0 = (1.0 / params[Acts::eBoundQOverP]) / Acts::UnitConstants::GeV;
            if (msgLevel(MSG::DEBUG)) {
              debug() << "track predicted p = " << p0 << " GeV" << endmsg;
            }
            if (std::abs(p0) > 500) {
              if (msgLevel(MSG::DEBUG)) {
                debug() << "skipping" << endmsg;
              }
              return;
            }

            eic::BasicParticle p{
                {-1, 0},
                eicd::sphericalToVector(1.0 / std::abs(params[Acts::eBoundQOverP]), 
                                               params[Acts::eBoundTheta],
                                               params[Acts::eBoundPhi]),
                {0., 0., 0.},                                                        // vectex 3-vector
                0.,                                                                  // time
                0,                                                                   // PDG particle code
                0,                                                                   // status
                static_cast<int16_t>(std::copysign(1., params[Acts::eBoundQOverP])), // charge
                1.                                                                   // weight
            };                                                                       // charge
            rec_parts->push_back(p);
          });
      }

      // set our IDs
      for (size_t i = 0; i < rec_parts->size(); ++i) {
        (*rec_parts)[i].ID({static_cast<int32_t>(i), 0});
      }
      
      return StatusCode::SUCCESS;
    }

  };
  DECLARE_COMPONENT(ParticlesFromTrackFit)

} // namespace Jug::Reco
