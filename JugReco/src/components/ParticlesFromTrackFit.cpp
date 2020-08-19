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
          debug() << trackstate.hasPredicted() << endmsg;
          debug() << trackstate.predicted() << endmsg;
          auto params = trackstate.predicted() ;//<< endmsg;
          debug() << 1.0/params[Acts::eQOP] << " GeV" << endmsg;

          eic::Particle p( {params[Acts::ePHI],params[Acts::eTHETA],1.0/params[Acts::eQOP],0.105}, {0.0,0.0,0.0,params[Acts::eT]},
                          (long long)13*params[Acts::eQOP]/std::abs(params[Acts::eQOP]), 0);
          debug() << p << endmsg;
          rec_parts->push_back(p);
        });

        //auto pos = m_geoSvc->cellIDPositionConverter()->position(ahit.cellID());
        //auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(ahit.cellID());
        //debug() << " dim size : " <<  std::size(dim) << endmsg;
        //for(const auto& s : dim ) {
        //  debug() << "a size : " <<  s << endmsg;
        //}
        ////std::array<double,3> posarr; pos.GetCoordinates(posarr);
        ////std::array<double,3> dimarr; dim.GetCoordinates(posarr);
        ////eic::TrackerHit hit;
        //eic::TrackerHit hit((long long)ahit.cellID(), (long long)ahit.cellID(), (long long)ahit.time(),
        //                    (float)ahit.charge() / 10000.0, (float)0.0, {{pos.x(), pos.y(),pos.z()}},{{dim[0],dim[1],0.0}});
        //rec_hits->push_back(hit);
      }
      return StatusCode::SUCCESS;
    }

  };
  DECLARE_COMPONENT(ParticlesFromTrackFit)

  } // namespace Examples
} // namespace Gaudi

