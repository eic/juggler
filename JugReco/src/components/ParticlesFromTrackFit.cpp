#include <algorithm>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Property.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
//#include "GaudiExamples/MyTrack.h"
#include "eicd/ParticleCollection.h"
#include "eicd/TrackerHitCollection.h"
#include "JugReco/SourceLinks.h"
#include "JugReco/Track.hpp"

namespace Jug {
  namespace Reco {
  
    /** Ultra-fast silicon detector digitization.
     *
     */
   class ParticlesFromTrackFit : public GaudiAlgorithm {
   public:
    //DataHandle<eic::RawTrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoryContainer>      m_inputTrajectories{"inputTrajectories", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer, this};

   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ParticlesFromTrackFit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputTrajectories", m_inputTrajectories,"");
          declareProperty("outputParticles", m_outputParticles, "");
        }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const TrajectoryContainer* trajectories = m_inputTrajectories.get();
      // create output collections
      auto rec_parts = m_outputParticles.createAndPut();

      for(const auto& traj : *trajectories) {
        //traj.trajectory().first
        auto tsize = traj.trajectory().first.size();
        debug() << "# fitted parameters : " << tsize << endmsg;
        if(tsize == 0 ) continue;
        traj.trajectory().second.visitBackwards(tsize-1, [&](auto&& trackstate) {
          //debug() << trackstate.hasPredicted() << endmsg;
          //debug() << trackstate.predicted() << endmsg;
          auto params = trackstate.predicted() ;//<< endmsg;

          double p0 = (1.0 / params[Acts::eBoundQOverP]) / Acts::UnitConstants::GeV;
          debug() << "track predicted p = " << p0 << " GeV" << endmsg;
          if ( std::abs(p0) > 500) {
            debug() << "skipping" << endmsg;
            return;
          }

          eic::Particle p({params[Acts::eBoundPhi], params[Acts::eBoundTheta], 1.0 / std::abs(params[Acts::eBoundQOverP]), 0.000511},
                          {0.0, 0.0, 0.0, params[Acts::eBoundTime]},
                          (long long)11 * params[Acts::eBoundQOverP] / std::abs(params[Acts::eBoundQOverP]), 0);
          //debug() << p << endmsg;
          rec_parts->push_back(p);
        });

      }
      return StatusCode::SUCCESS;
    }

  };
  DECLARE_COMPONENT(ParticlesFromTrackFit)

  } // namespace Examples
} // namespace Gaudi

