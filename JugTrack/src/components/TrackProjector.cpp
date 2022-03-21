#include <algorithm>

// FIXME needs renaming (TrackProjector) and updating to fix and remove hardcoded numbers

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
#include "eicd/TrajectoryCollection.h"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Trajectories.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

#include "eicd/VectorPolar.h"

#include <cmath>

namespace Jug::Reco {

  /** Extrac the particles form fit trajectories.
   *
   * \ingroup tracking
   */
   class TrajectoryFromTrackFit : public GaudiAlgorithm, AlgorithmIDMixin<int32_t> {
   public:
    DataHandle<TrajectoriesContainer>     m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<eicd::TrajectoryCollection> m_outputTrajectory{"outputTrajectoryParameters", Gaudi::DataHandle::Writer, this};

   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    TrajectoryFromTrackFit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
        , AlgorithmIDMixin(name, info()) {
          declareProperty("inputTrajectories", m_inputTrajectories,"");
          declareProperty("outputTrajectoryParameters", m_outputTrajectory, "ACTS Trajectory Parameters");
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
      auto traj_pars = m_outputTrajectory.createAndPut();

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
        debug() << "# of elements in trackTips " <<trackTips.size() << endmsg;

        if (trackTips.empty()) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "Empty multiTrajectory." << endmsg;
          }
          continue;
        }

        auto& trackTip = trackTips.front();

        // Collect the trajectory summary info
        auto trajState       = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
        int  m_nMeasurements = trajState.nMeasurements;
        int  m_nStates       = trajState.nStates;
        int  m_nCalibrated   = 0;
        if (msgLevel(MSG::DEBUG)) {
    		  debug() << "n measurement in trajectory " << m_nMeasurements << endmsg;
    		  debug() << "n state in trajectory " << m_nStates << endmsg;
        }

        // get path length at last silicon layer
        float pathlength_at_reflayer = -9999.;
        mj.visitBackwards(trackTip, [&](auto&& trackstate) {
          // debug() << trackstate.hasPredicted() << endmsg;
          // debug() << trackstate.predicted() << endmsg;
          auto params = trackstate.predicted(); //<< endmsg;
          auto pathlength = trackstate.pathLength();
          auto geoID = trackstate.referenceSurface().geometryId();
          auto volume = geoID.volume();
          auto layer = geoID.layer();
          if (trackstate.hasCalibrated())
          {
            m_nCalibrated++;
          } 

          if (pathlength_at_reflayer<0 && volume==29 && layer==4 && trackstate.hasCalibrated()) pathlength_at_reflayer = pathlength; // 2nd outter barrel
          if (pathlength_at_reflayer<0 && pathlength>1700) pathlength_at_reflayer = pathlength; // endcap GEM layer

          if (msgLevel(MSG::DEBUG)) {
            debug() << "******************************" << endmsg;
            debug() << "predicted variables: \n" << trackstate.predicted() << endmsg;
            debug() << "geoID = " << geoID << endmsg;
            debug() << "volume = " << volume << ", layer = " << layer << endmsg;
            debug() << "pathlength = " << pathlength << endmsg;
            debug() << "hasCalibrated = " << trackstate.hasCalibrated() << endmsg;
            debug() << "******************************" << endmsg;
          }
        });

        if (msgLevel(MSG::DEBUG)) {
          debug() << "n calibrated state in trajectory " << m_nCalibrated << endmsg;
        }

        // Get the fitted track parameter
        bool m_hasFittedParams = false;
        int ID = 0;
        float charge = -9999.;
        float tof = -9999.;
        if (traj.hasTrackParameters(trackTip)) {
          m_hasFittedParams      = true;
          const auto& boundParam = traj.trackParameters(trackTip);
          const auto& parameter  = boundParam.parameters();
          const auto& covariance = *boundParam.covariance();
          charge = boundParam.charge();

          double p = 1.0 / parameter[Acts::eBoundQOverP];
          double p_phi = parameter[Acts::eBoundPhi];
          double p_theta = parameter[Acts::eBoundTheta];
          eicd::VectorXYZT mom(p*sin(p_theta)*cos(p_phi), p*sin(p_theta)*sin(p_phi), p*cos(p_theta), 0);

          if (msgLevel(MSG::DEBUG)) {
            debug() << " chi2 = " << trajState.chi2Sum << endmsg;
            debug() << " pathlength at reference layer = " << pathlength_at_reflayer << "mm" << endmsg;
          }
            
          eicd::Trajectory traj_par{
            {ID++, algorithmID()},
            {0, 0}, // proto track ID
            {0, 0}, // track param ID
            mom, // vectex 4-vector
            pathlength_at_reflayer,
            charge,
            tof}; 

          traj_pars->push_back(traj_par);
        }
      }
      
      return StatusCode::SUCCESS;
    }

  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(TrajectoryFromTrackFit)

} // namespace Jug::Reco
