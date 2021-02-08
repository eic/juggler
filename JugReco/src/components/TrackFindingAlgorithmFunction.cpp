#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

//#include "ActsExamples/Plugins/BField/ScalableBField.hpp"
//#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include "TrackFindingAlgorithm.h"

#include "JugReco/BField.h"

#include <random>
#include <stdexcept>

namespace {
  template <typename track_finder_t>
  struct TrackFinderFunctionImpl {
    track_finder_t trackFinder;

    TrackFinderFunctionImpl(track_finder_t&& f) : trackFinder(std::move(f)) {}

    Jug::Reco::TrackFindingAlgorithm::TrackFinderResult
    operator()(const Jug::IndexSourceLinkContainer&                  sourcelinks,
               const Jug::TrackParametersContainer&                  initialParameters,
               const Jug::Reco::TrackFindingAlgorithm::TrackFinderOptions& options) const
    {
      return trackFinder.findTracks(sourcelinks, initialParameters, options);
    };
  };
} // namespace

//ActsExamples::TrackFindingAlgorithm::TrackFinderFunction ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction(
//    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry, Options::BFieldVariant magneticField)
//{
//  using Updater  = Acts::GainMatrixUpdater;
//  using Smoother = Acts::GainMatrixSmoother;
//
//  // unpack the magnetic field variant and instantiate the corresponding track
//  // finder.
//  return std::visit(
//      [trackingGeometry](auto&& inputField) -> TrackFinderFunction {
//        // each entry in the variant is already a shared_ptr
//        // need ::element_type to get the real magnetic field type
//        using InputMagneticField =
//            typename std::decay_t<decltype(inputField)>::element_type;
//        using MagneticField = Acts::SharedBField<InputMagneticField>;
//        using Stepper = Acts::EigenStepper<MagneticField>;
//        using Navigator = Acts::Navigator;
//        using Propagator = Acts::Propagator<Stepper, Navigator>;
//        using CKF =
//            Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother>;
//
//        // construct all components for the track finder
//        MagneticField field(std::move(inputField));
//        Stepper stepper(std::move(field));
//        Navigator navigator(trackingGeometry);
//        navigator.resolvePassive = false;
//        navigator.resolveMaterial = true;
//        navigator.resolveSensitive = true;
//        Propagator propagator(std::move(stepper), std::move(navigator));
//        CKF trackFinder(std::move(propagator));
//
//        // build the track finder functions. owns the track finder object.
//        return TrackFinderFunctionImpl<CKF>(std::move(trackFinder));
//      },
//      std::move(magneticField));
//}

namespace Jug::Reco {

  TrackFindingAlgorithm::TrackFinderFunction
  TrackFindingAlgorithm::makeTrackFinderFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                                 BFieldVariant                                 magneticField)
  {
    using Updater  = Acts::GainMatrixUpdater;
    using Smoother = Acts::GainMatrixSmoother;

    // unpack the magnetic field variant and instantiate the corresponding track
    // finder.
    return std::visit(
        [trackingGeometry](auto&& inputField) -> TrackFinderFunction {
          // each entry in the variant is already a shared_ptr
          // need ::element_type to get the real magnetic field type
          using InputMagneticField = typename std::decay_t<decltype(inputField)>::element_type;
          using MagneticField      = Acts::SharedBField<InputMagneticField>;
          using Stepper            = Acts::EigenStepper<MagneticField>;
          using Navigator          = Acts::Navigator;
          using Propagator         = Acts::Propagator<Stepper, Navigator>;
          using CKF                = Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother>;

          //std::cout << " finding ...\n";
          // construct all components for the track finder
          MagneticField field(std::move(inputField));
          Stepper       stepper(std::move(field));
          Navigator     navigator(trackingGeometry);
          navigator.resolvePassive   = false;
          navigator.resolveMaterial  = true;
          navigator.resolveSensitive = true;
          //std::cout << " propagator\n";
          Propagator propagator(std::move(stepper), std::move(navigator));
          CKF        trackFinder(std::move(propagator));

          // build the track finder functions. owns the track finder object.
          return TrackFinderFunctionImpl<CKF>(std::move(trackFinder));
        },
        std::move(magneticField));
  }

} // namespace Jug::Reco

