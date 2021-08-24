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
  using Fitter           = Acts::KalmanFitter<Propagator, Updater, Smoother>;
  using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
  using DirectFitter     = Acts::KalmanFitter<DirectPropagator, Updater, Smoother>;


  /** Fitter implmentation .
   *
   * \ingroup tracking
   */
  template <typename track_fitter_t>
  struct TrackFitterFunctionImpl {
    track_fitter_t trackFitter;

    TrackFitterFunctionImpl(track_fitter_t&& f) : trackFitter(std::move(f)) {}

    Jug::Reco::TrackFittingAlgorithm::FitterResult
    operator()(const std::vector<Jug::IndexSourceLink>&                    sourceLinks,
               const Jug::TrackParameters&                                 initialParameters,
               const Jug::Reco::TrackFittingAlgorithm::TrackFitterOptions& options) const
    {
      return trackFitter.fit(sourceLinks, initialParameters, options);
    }

    //TrackFittingAlgorithm::TrackFitterResult
    //operator()(const std::vector<SourceLink>&              sourceLinks,
    //           const TrackParameters&                           initialParameters,
    //           const TrackFittingAlgorithm::TrackFitterOptions& options) const
    //{
    //  return trackFitter.fit(sourceLinks, initialParameters, options);
    //};
  };

}  // namespace

namespace Jug::Reco {

  //using InputMagneticField = typename std::decay_t<decltype(inputField)>::element_type;
  //using MagneticField      = Acts::SharedBField<InputMagneticField>;
  //using Stepper            = Acts::EigenStepper<MagneticField>;
  //using Navigator          = Acts::Navigator;
  //using Propagator         = Acts::Propagator<Stepper, Navigator>;
  //using Fitter             = Acts::KalmanFitter<Propagator, Updater, Smoother>;

  using Updater          = Acts::GainMatrixUpdater;
  using Smoother         = Acts::GainMatrixSmoother;
  using Stepper          = Acts::EigenStepper<>;
  using Propagator       = Acts::Propagator<Stepper, Acts::Navigator>;
  using Fitter           = Acts::KalmanFitter<Propagator, Updater, Smoother>;
  using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
  using DirectFitter     = Acts::KalmanFitter<DirectPropagator, Updater, Smoother>;

  TrackFittingAlgorithm::FitterFunction TrackFittingAlgorithm::makeTrackFittingFunction(
      std::shared_ptr<const Acts::TrackingGeometry>      trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField)

  {
    Stepper                 stepper(std::move(magneticField));
    //Acts::Navigator::Config cfg{trackingGeometry};
    //cfg.resolvePassive   = false;
    //cfg.resolveMaterial  = true;
    //cfg.resolveSensitive = true;
    Acts::Navigator navigator(trackingGeometry);
    navigator.resolvePassive   = false;
    navigator.resolveMaterial  = true;
    navigator.resolveSensitive = true;
    Propagator      propagator(std::move(stepper), std::move(navigator));
    Fitter          trackFitter(std::move(propagator));

    // build the fitter functions. owns the fitter object.
    return TrackFitterFunctionImpl(std::move(trackFitter));
  }
  //  using Updater  = Acts::GainMatrixUpdater;
  //  using Smoother = Acts::GainMatrixSmoother;

  //  // unpack the magnetic field variant and instantiate the corresponding fitter.
  //  return std::visit(
  //      [trackingGeometry](auto&& inputField) -> FitterFunction {
  //        // each entry in the variant is already a shared_ptr
  //        // need ::element_type to get the real magnetic field type
  //        // construct all components for the fitter
  //        MagneticField field(std::move(inputField));
  //        Stepper       stepper(std::move(field));
  //        Navigator     navigator(trackingGeometry);
  //        navigator.resolvePassive   = false;
  //        navigator.resolveMaterial  = true;
  //        navigator.resolveSensitive = true;
  //        Propagator propagator(std::move(stepper), std::move(navigator));
  //        Fitter     trackFitter(std::move(propagator));

  //        // build the fitter functions. owns the fitter object.
  //        return TrackFitterFunctionImpl<Fitter>(std::move(trackFitter));
  //      },
  //      std::move(magneticField));
  //}
} // namespace Jug::Reco
