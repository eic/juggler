// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Whitney Armstrong
#include <algorithm>

// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

// Event Model related classes
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Trajectories.hpp"
#include "eicd/BasicParticleCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/TrackerHitCollection.h"

#include "Acts/Utilities/Helpers.hpp"

#include "eicd/VectorPolar.h"

#include <cmath>

namespace Jug::Track {

/** Read and convert the Acts trajectory information to create EICD trajectory and track
 *  parameter structures
 *
 * \ingroup tracking
 */
class ActsTrajectoryReader : public GaudiAlgorithm {
private:
  DataHandle<eicd::TrackerHitCollection> m_inputHits{"inputHits", Gaudi::DataHandle::Reader, this};
  DataHandle<TrajectoriesContainer> m_inputTrajectories{"inputActsTrajectories", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::TrajectoryCollection> m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};
  DataHandle<eicd::TrackParametersCollection> m_outputParameters{"outputParameters", Gaudi::DataHandle::Writer, this};

public:
  //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
  ActsTrajectoryReader(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHits", m_inputHits, "");
    declareProperty("inputActsTrajectories", m_inputTrajectories, "");
    declareProperty("outputTrajectories", m_outputTrajectories, "");
    declareProperty("outputParameters", m_outputParameters, "Acts Track Parameters");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input
    const TrajectoriesContainer& acts_traj = *(m_inputTrajectories.get());
    const eicd::TrackerHitCollection& hits  = *(m_inputHits.get());
    // create output collections
    auto out_traj = m_outputTrajectories.createAndPut();
    auto out_pars = m_outputParameters.createAndPut();

    if (msgLevel(MSG::DEBUG)) {
      debug() << std::size(*trajectories) << " trajectories " << endmsg;
    }

    // Loop over the trajectories
    for (const auto& traj : acts_traj) {

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
      auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
      // int  m_nMeasurements = trajState.nMeasurements;
      // int  m_nStates       = trajState.nStates;

      // Get the fitted track parameter
      //
      bool hasFittedParams = false;
      if (traj.hasTrackParameters(trackTip)) {
        hasFittedParams        = true;
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

        const std::array<float, 21> covMatrix{static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)),
                                              static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)),
                                              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
                                              static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundPhi)),
                                              static_cast<float>(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)),
                                              static_cast<float>(covariance(Acts::eBoundTime, Acts::eBoundTime)),
                                              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc1)),
                                              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundTheta)),
                                              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundPhi)),
                                              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundQOverP)),
                                              static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundTime)),
                                              static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundTheta)),
                                              static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundPhi)),
                                              static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundQOverP)),
                                              static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundTime)),
                                              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundPhi)),
                                              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundQOverP)),
                                              static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundTime)),
                                              static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundQOverP)),
                                              static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundTime)),
                                              static_cast<float>(covariance(Acts::eBoundQOverP, Acts::eBoundTime))};

        eicd::TrackParameters pars{0, // type: track head --> 0
                                  {parameter[Acts::eBoundLoc0], parameter[Acts::eBoundLoc1]},
                                  static_cast<float>(parameter[Acts::eBoundTheta]),
                                  static_cast<float>(parameter[Acts::eBoundPhi]),
                                  static_cast<float>(parameter[Acts::eBoundQOverP]),
                                  static_cast<float>(parameter[Acts::eBoundTime]),
                                  covMatrix};
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

        eicd::BasicParticle p{
            {-1, 0},
            eicd::VectorPolar( // 3-momentum vector
                {1.0 / std::abs(params[Acts::eBoundQOverP]), params[Acts::eBoundTheta], params[Acts::eBoundPhi]}),
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
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ActsTrajectoryReader)

} // namespace Jug::Track
