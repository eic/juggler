// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten

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

#include "TrackFindingAlgorithm.h"

#include "JugBase/BField/DD4hepBField.h"


#include <random>
#include <stdexcept>

namespace {
  using Updater  = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  using Stepper    = Acts::EigenStepper<>;
  using Navigator  = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using CKF        = Acts::CombinatorialKalmanFilter<Propagator>;

  /** Finder implmentation .
   *
   * \ingroup track
   */
  struct TrackFinderFunctionImpl
    : public Jug::Reco::TrackFindingAlgorithm::TrackFinderFunction {
    CKF trackFinder;

    TrackFinderFunctionImpl(CKF&& f) : trackFinder(std::move(f)) {}

    Jug::Reco::TrackFindingAlgorithm::TrackFinderResult
    operator()(const Jug::IndexSourceLinkContainer&                        sourcelinks,
               const Jug::TrackParametersContainer&                        initialParameters,
               const Jug::Reco::TrackFindingAlgorithm::TrackFinderOptions& options)
               const override
    {
      return trackFinder.findTracks(sourcelinks, initialParameters, options);
    };
  };

} // namespace

namespace Jug::Reco {

  std::shared_ptr<TrackFindingAlgorithm::TrackFinderFunction>
  TrackFindingAlgorithm::makeTrackFinderFunction(
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
    return std::make_shared<TrackFinderFunctionImpl>(std::move(trackFinder));
  }

} // namespace Jug::Reco
