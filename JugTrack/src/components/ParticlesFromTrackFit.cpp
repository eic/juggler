// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck

#include <algorithm>

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

// Event Model related classes
#include "edm4eic/EDM4eicVersion.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/TrackerHitCollection.h"
#include "edm4eic/TrackParametersCollection.h"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include "edm4hep/utils/vector_utils.h"

#include <cmath>

namespace Jug::Reco {

  #if EDM4EIC_VERSION_MAJOR >= 5
    // This array relates the Acts and EDM4eic covariance matrices, including
    // the unit conversion to get from Acts units into EDM4eic units.
    //
    // Note: std::map is not constexpr, so we use a constexpr std::array
    // std::array initialization need double braces since arrays are aggregates
    // ref: https://en.cppreference.com/w/cpp/language/aggregate_initialization
    static constexpr std::array<std::pair<Acts::BoundIndices, double>, 6> edm4eic_indexed_units{{
      {Acts::eBoundLoc0, Acts::UnitConstants::mm},
      {Acts::eBoundLoc1, Acts::UnitConstants::mm},
      {Acts::eBoundTheta, 1.},
      {Acts::eBoundPhi, 1.},
      {Acts::eBoundQOverP, 1. / Acts::UnitConstants::GeV},
      {Acts::eBoundTime, Acts::UnitConstants::ns}
    }};
  #endif

  /** Extract the particles form fit trajectories.
   *
   * \ingroup tracking
   */
   class ParticlesFromTrackFit : public Gaudi::Algorithm {
   private:
    mutable DataHandle<ActsExamples::TrajectoriesContainer>     m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    mutable DataHandle<edm4eic::ReconstructedParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer, this};
    mutable DataHandle<edm4eic::TrackParametersCollection> m_outputTrackParameters{"outputTrackParameters", Gaudi::DataHandle::Writer, this};

   public:
    //  ill-formed: using Gaudi::Algorithm::GaudiAlgorithm;
    ParticlesFromTrackFit(const std::string& name, ISvcLocator* svcLoc)
        : Gaudi::Algorithm(name, svcLoc) {
          declareProperty("inputTrajectories", m_inputTrajectories,"");
          declareProperty("outputParticles", m_outputParticles, "");
          declareProperty("outputTrackParameters", m_outputTrackParameters, "Acts Track Parameters");
        }

    StatusCode initialize() override {
      if (Gaudi::Algorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute(const EventContext&) const override {
      // input collection
      const auto* const trajectories = m_inputTrajectories.get();
      // create output collections
      auto* rec_parts = m_outputParticles.createAndPut();
      auto* track_pars = m_outputTrackParameters.createAndPut();

      if (msgLevel(MSG::DEBUG)) {
        debug() << std::size(*trajectories) << " trajectories " << endmsg;
      }

      // Loop over the trajectories
        for (const auto& traj : *trajectories) {

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

          const auto& trackTip = trackTips.front();

          // Collect the trajectory summary info
          auto trajState       = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
          //int  m_nMeasurements = trajState.nMeasurements;
          //int  m_nStates       = trajState.nStates;

          // Get the fitted track parameter
          //
          if (traj.hasTrackParameters(trackTip)) {
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

            auto pars = track_pars->create();
            pars.setType(0); // type: track head --> 0
            pars.setLoc({
              static_cast<float>(parameter[Acts::eBoundLoc0]),
              static_cast<float>(parameter[Acts::eBoundLoc1])
            });
            pars.setTheta(static_cast<float>(parameter[Acts::eBoundTheta]));
            pars.setPhi(static_cast<float>(parameter[Acts::eBoundPhi]));
            pars.setQOverP(static_cast<float>(parameter[Acts::eBoundQOverP]));
            pars.setTime(static_cast<float>(parameter[Acts::eBoundTime]));
            #if EDM4EIC_VERSION_MAJOR >= 5
              edm4eic::Cov6f cov;
              for (size_t i = 0; const auto& [a, x] : edm4eic_indexed_units) {
                for (size_t j = 0; const auto& [b, y] : edm4eic_indexed_units) {
                  // FIXME why not pars.getCovariance()(i,j) = covariance(a,b) / x / y;
                  cov(i,j) = covariance(a,b) / x / y;
                }
              }
              pars.setCovariance(cov);
            #else
              pars.setCharge(static_cast<float>(boundParam.charge()));
              pars.setLocError({
                    static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)),
                    static_cast<float>(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)),
                    static_cast<float>(covariance(Acts::eBoundLoc0, Acts::eBoundLoc1))
                });
              pars.setMomentumError({
                    static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
                    static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundPhi)),
                    static_cast<float>(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)),
                    static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundPhi)),
                    static_cast<float>(covariance(Acts::eBoundTheta, Acts::eBoundQOverP)),
                    static_cast<float>(covariance(Acts::eBoundPhi, Acts::eBoundQOverP))
                });
              pars.setTimeError(sqrt(static_cast<float>(covariance(Acts::eBoundTime, Acts::eBoundTime))));
            #endif
          }

          auto tsize = trackTips.size();
          if (msgLevel(MSG::DEBUG)) {
            debug() << "# fitted parameters : " << tsize << endmsg;
          }
          if (tsize == 0) {
            continue;
          }

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

            auto rec_part = rec_parts->create();
            rec_part.setMomentum(
              edm4hep::utils::sphericalToVector(
                1.0 / std::abs(params[Acts::eBoundQOverP]),
                params[Acts::eBoundTheta],
                params[Acts::eBoundPhi])
            );
            rec_part.setCharge(static_cast<int16_t>(std::copysign(1., params[Acts::eBoundQOverP])));
          });
      }

      return StatusCode::SUCCESS;
    }

  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(ParticlesFromTrackFit)

} // namespace Jug::Reco
