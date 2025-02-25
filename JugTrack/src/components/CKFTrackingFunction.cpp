// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

#if 0
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#endif

#include "CKFTracking.h"

#include "JugTrack/DD4hepBField.h"


#include <random>
#include <stdexcept>

namespace {
  using Updater  = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  using Stepper    = Acts::EigenStepper<>;
  using Navigator  = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
#if Acts_VERSION_MAJOR >= 36
  using CKF =
      Acts::CombinatorialKalmanFilter<Propagator, ActsExamples::TrackContainer>;
#else
  using CKF =
      Acts::CombinatorialKalmanFilter<Propagator, Acts::VectorMultiTrajectory>;

  using TrackContainer =
      Acts::TrackContainer<Acts::VectorTrackContainer,
                           Acts::VectorMultiTrajectory, std::shared_ptr>;
#endif

  /** Finder implmentation .
   *
   * \ingroup track
   */
  struct CKFTrackingFunctionImpl
    : public Jug::Reco::CKFTracking::CKFTrackingFunction {
    CKF trackFinder;

    CKFTrackingFunctionImpl(CKF&& f) : trackFinder(std::move(f)) {}

    Jug::Reco::CKFTracking::TrackFinderResult
    operator()(const ActsExamples::TrackParameters& initialParameters,
               const Jug::Reco::CKFTracking::TrackFinderOptions& options,
#if Acts_VERSION_MAJOR >= 36
               ActsExamples::TrackContainer& tracks) const override {
#else
               TrackContainer& tracks) const override {
#endif
      return trackFinder.findTracks(initialParameters, options, tracks);
    };
  };

} // namespace

namespace Jug::Reco {

  std::shared_ptr<CKFTracking::CKFTrackingFunction>
  CKFTracking::makeCKFTrackingFunction(
      std::shared_ptr<const Acts::TrackingGeometry>      trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField)
  {
    Stepper   stepper(std::move(magneticField));
    Navigator::Config cfg{trackingGeometry};
    cfg.resolvePassive   = false;
    cfg.resolveMaterial  = true;
    cfg.resolveSensitive = true;
    Navigator navigator(cfg);

    Propagator propagator(std::move(stepper), std::move(navigator));
    CKF        trackFinder(std::move(propagator));

    // build the track finder functions. onws the track finder object.
    return std::make_shared<CKFTrackingFunctionImpl>(std::move(trackFinder));
  }

} // namespace Jug::Reco
