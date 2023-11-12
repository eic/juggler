// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "TrackFittingAlgorithm.h"

namespace {

  using Updater          = Acts::GainMatrixUpdater;
  using Smoother         = Acts::GainMatrixSmoother;
  using Stepper          = Acts::EigenStepper<>;
  using Propagator       = Acts::Propagator<Stepper, Acts::Navigator>;
  using Fitter           = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;
  using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
  using DirectFitter     = Acts::KalmanFitter<DirectPropagator, Acts::VectorMultiTrajectory>;

  /** Fitter implmentation .
   *
   * \ingroup tracking
   */
  template <typename track_fitter_t>
  struct TrackFitterFunctionImpl {
    track_fitter_t trackFitter;

    TrackFitterFunctionImpl(track_fitter_t&& f) : trackFitter(std::move(f)) {}

    Jug::Reco::TrackFittingAlgorithm::FitterResult
    operator()(const std::vector<ActsExamples::IndexSourceLink>&                    sourceLinks,
               const ActsExamples::TrackParameters&                                 initialParameters,
               const Jug::Reco::TrackFittingAlgorithm::TrackFitterOptions& options) const
    {
      return trackFitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters, options);
    }

  };

}  // namespace

namespace Jug::Reco {

  using Updater          = Acts::GainMatrixUpdater;
  using Smoother         = Acts::GainMatrixSmoother;
  using Stepper          = Acts::EigenStepper<>;
  using Propagator       = Acts::Propagator<Stepper, Acts::Navigator>;
  using Fitter           = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;
  using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
  using DirectFitter     = Acts::KalmanFitter<DirectPropagator, Acts::VectorMultiTrajectory>;

  TrackFittingAlgorithm::FitterFunction TrackFittingAlgorithm::makeTrackFittingFunction(
      std::shared_ptr<const Acts::TrackingGeometry>      trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField)

  {
    Stepper                 stepper(std::move(magneticField));
    Acts::Navigator::Config cfg{trackingGeometry};
    cfg.resolvePassive   = false;
    cfg.resolveMaterial  = true;
    cfg.resolveSensitive = true;
    Acts::Navigator navigator(cfg);
    Propagator      propagator(std::move(stepper), std::move(navigator));
    Fitter          trackFitter(std::move(propagator));

    // build the fitter functions. owns the fitter object.
    return TrackFitterFunctionImpl(std::move(trackFitter));
  }

} // namespace Jug::Reco
