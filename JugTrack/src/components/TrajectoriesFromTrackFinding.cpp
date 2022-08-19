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

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

// Event Model related classes
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

#include <cmath>

namespace Jug::Reco {

  /** Extrac the particles form fit trajectories.
   *
   * \ingroup tracking
   */
   class TrajectoriesFromTrackFinding : public GaudiAlgorithm {
   public:
      DataHandle<TrajectoriesContainer>     m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::TrajectoryCollection> m_outputTrajectory{"outputTrajectoryParameters", Gaudi::DataHandle::Writer, this};

   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
   TrajectoriesFromTrackFinding(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
        {
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

      // ActsExamples::ProcessCode ActsExamples::RootTrajectoryStatesWriter::writeT(
    // const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
    // using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
    // using HitSimHitsMap = IndexMultimap<Index>;

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
          if (trackTips.empty()) {
            if (msgLevel(MSG::DEBUG)) {
              debug() << "Empty multiTrajectory." << endmsg;
            }
            continue;
          }

          auto& trackTip = trackTips.front();

          // Collect the trajectory summary info
          auto  trajState       = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
          unsigned int  m_nStates       = trajState.nStates;
          unsigned int  m_nMeasurements = trajState.nMeasurements;
          unsigned int  m_nOutliers     = trajState.nOutliers;
          unsigned int  m_nHoles        = trajState.nHoles;
          float         m_chi2Sum       = trajState.chi2Sum;
          unsigned int  m_ndf           = trajState.NDF;    
          unsigned int  m_nSharedHits   = trajState.nSharedHits;  

          unsigned int hasFit=0;
          if (traj.hasTrackParameters(trackTip))
            hasFit=1;

          if (msgLevel(MSG::DEBUG)) {

              debug() << "hasFit                       = " << hasFit<<endmsg;
              debug() << "states                       = " << m_nStates<<endmsg;
              debug() << "Measurements                 = " << m_nMeasurements<<endmsg;
              debug() << "Holes, outliers, shared hits = " << m_nHoles<<"/"<<m_nOutliers<<"/"<<m_nSharedHits<<endmsg;
              debug() << "Total chi2/ndf               = " << m_chi2Sum <<"/"<< m_ndf<< endmsg;
            }

          eic::Trajectory traj_par{hasFit,
            m_nStates,
            m_nMeasurements,
            m_nOutliers, 
            m_nHoles,
            m_chi2Sum, 
            m_ndf, 
            m_nSharedHits
            // trajState.measurementChi2,
            // trajState.outlierChi2, 
            // 0,0
          };

          for(unsigned int ii=0;ii<m_nMeasurements;ii++){
            float cc = trajState.measurementChi2[ii];
            traj_par.addmeasurementChi2(cc);
            if (msgLevel(MSG::DEBUG)) 
              debug() << "measurementChi2/vol/layer " << ii<<":"<<cc<<"/"<<trajState.measurementVolume[ii]<<"/"<<trajState.measurementLayer[ii]<<"/"<<endmsg;
          }

          for(unsigned int ii=0;ii<m_nOutliers;ii++){
            float cc = trajState.outlierChi2[ii];
            traj_par.addoutlierChi2(cc);
            if (msgLevel(MSG::DEBUG)) 
              debug() << "outlierChi2/vol/layer " << ii<<":"<<cc<<"/"<<trajState.outlierVolume[ii]<<"/"<<trajState.outlierLayer[ii]<<"/"<<endmsg;
          }

          traj_pars->push_back(traj_par);

//
// corresponding eicd structure
//
    // Description: "Raw trajectory from the tracking algorithm"
    // Author: "S. Joosten, S. Li"
    // Members:
    //   - uint32_t          type              // 0 (does not have good track fit), 1 (has good track fit)
    //   - uint32_t          nStates           // Number of tracking steps
    //   - uint32_t          nMeasurements     // Number of hits used 
    //   - uint32_t          nOutliers         // Number of hits not considered 
    //   - uint32_t          nHoles            // Number of missing hits
    //   - float             chi2              // Total chi2
    //   - uint32_t          ndf               // Number of degrees of freedom
    //   - uint32_t          nSharedHits       // Number of shared hits with other trajectories
    // VectorMembers:
    //   - float             measurementChi2   // Chi2 for each of the measurements
    //   - float             outlierChi2       // Chi2 for each of the outliers
    // OneToManyRelations:
    //   - eic::TrackerHit   measurementHits   // Measurement hits used in this trajectory
    //   - eic::TrackerHit   outlierHits       // Outlier hits not used in this trajectory

//
//// to do: get hits info (e.g. position)
//
      // // Get the trackStates on the trajectory
      // m_nParams = {0, 0, 0};
      // mj.visitBackwards(trackTip, [&](const auto& state) {
      //   // we only fill the track states with non-outlier measurement
      //   auto typeFlags = state.typeFlags();
      //   if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      //     return true;
      //   }

      //   const auto& surface = state.referenceSurface();

      //   // get the truth hits corresponding to this trackState
      //   // Use average truth in the case of multiple contributing sim hits
      //   const auto& sl =
      //       static_cast<const IndexSourceLink&>(state.uncalibrated());
      //   const auto hitIdx = sl.index();
      //   auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
      //   auto [truthLocal, truthPos4, truthUnitDir] =
      //       averageSimHits(ctx.geoContext, surface, simHits, indices);
      //   // momemtum averaging makes even less sense than averaging position and
      //   // direction. use the first momentum or set q/p to zero
      //   float truthQOP = 0.0f;
      //   if (not indices.empty()) {
      //     // we assume that the indices are within valid ranges so we do not
      //     // need to check their validity again.
      //     const auto simHitIdx0 = indices.begin()->second;
      //     const auto& simHit0 = *simHits.nth(simHitIdx0);
      //     const auto p =
      //         simHit0.momentum4Before().template segment<3>(Acts::eMom0).norm();
      //     truthQOP = truthQ / p;
      //   }

      //   // fill the truth hit info
      //   m_t_x.push_back(truthPos4[Acts::ePos0]);
      //   m_t_y.push_back(truthPos4[Acts::ePos1]);
      //   m_t_z.push_back(truthPos4[Acts::ePos2]);
      //   m_t_r.push_back(perp(truthPos4.template segment<3>(Acts::ePos0)));
      //   m_t_dx.push_back(truthUnitDir[Acts::eMom0]);
      //   m_t_dy.push_back(truthUnitDir[Acts::eMom1]);
      //   m_t_dz.push_back(truthUnitDir[Acts::eMom2]);

      //   // get the truth track parameter at this track State
      //   float truthLOC0 = truthLocal[Acts::ePos0];
      //   float truthLOC1 = truthLocal[Acts::ePos1];
      //   float truthTIME = truthPos4[Acts::eTime];
      //   float truthPHI = phi(truthUnitDir);
      //   float truthTHETA = theta(truthUnitDir);

      //   // fill the truth track parameter at this track State
      //   m_t_eLOC0.push_back(truthLOC0);
      //   m_t_eLOC1.push_back(truthLOC1);
      //   m_t_ePHI.push_back(truthPHI);
      //   m_t_eTHETA.push_back(truthTHETA);
      //   m_t_eQOP.push_back(truthQOP);
      //   m_t_eT.push_back(truthTIME);

      //   // get the geometry ID
      //   auto geoID = surface.geometryId();
      //   m_volumeID.push_back(geoID.volume());
      //   m_layerID.push_back(geoID.layer());
      //   m_moduleID.push_back(geoID.sensitive());


      }
      return StatusCode::SUCCESS;
    }

  };
  DECLARE_COMPONENT(TrajectoriesFromTrackFinding)

} // namespace Jug::Reco
