#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

#if 0
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#endif

#include "TrackFindingAlgorithm.h"

#include "JugTrack/BField.h"

#include <random>
#include <stdexcept>

namespace {
  using Updater  = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  using Stepper    = Acts::EigenStepper<>;
  using Navigator  = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using CKF        = Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother>;

  struct TrackFinderFunctionImpl {
    CKF trackFinder;

    TrackFinderFunctionImpl(CKF&& f) : trackFinder(std::move(f)) {}

    Jug::Reco::TrackFindingAlgorithm::TrackFinderResult
    operator()(const Jug::IndexSourceLinkContainer&                        sourcelinks,
               const Jug::TrackParametersContainer&                        initialParameters,
               const Jug::Reco::TrackFindingAlgorithm::TrackFinderOptions& options) const
    {
      return trackFinder.findTracks(sourcelinks, initialParameters, options);
    };
  };

} // namespace

namespace Jug::Reco {

  TrackFindingAlgorithm::TrackFinderFunction TrackFindingAlgorithm::makeTrackFinderFunction(
      std::shared_ptr<const Acts::TrackingGeometry>      trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField)
  {
    Stepper   stepper(std::move(magneticField));
    Navigator navigator(trackingGeometry);
    navigator.resolvePassive   = false;
    navigator.resolveMaterial  = true;
    navigator.resolveSensitive = true;

    Propagator propagator(std::move(stepper), std::move(navigator));
    CKF        trackFinder(std::move(propagator));

    // build the track finder functions. onws the track finder object.
    return TrackFinderFunctionImpl(std::move(trackFinder));
  }

} // namespace Jug::Reco
