#include <cmath>
// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugReco/Track.hpp"
#include "Acts/Utilities/Units.hpp"

#include "eicd/TrackerHitCollection.h"
#include "eicd/ClusterCollection.h"


  ///// (Reconstructed) track parameters e.g. close to the vertex.
  //using TrackParameters = Acts::CurvilinearTrackParameters;

  ///// Container of reconstructed track states for multiple tracks.
  //using TrackParametersContainer = std::vector<TrackParameters>;

  ///// MultiTrajectory definition
  //using Trajectory = Acts::MultiTrajectory<SourceLink>;

  ///// Container for the truth fitting/finding track(s)
  //using TrajectoryContainer = std::vector<SimMultiTrajectory>;

namespace Jug::Reco {

  /** Initial Track parameters from clusters and vertex hits.
   *
   *
   * The momentum of the initial track is estimated from the cluster  energy and 
   * the direction is set using the vertex hits.
   *
   */
  class TrackParamVertexClusterInit : public GaudiAlgorithm {
  public:
    using Clusters =  eic::ClusterCollection;
    using VertexHits =  eic::TrackerHitCollection;

    DataHandle<VertexHits>               m_inputVertexHits{"inputVertexHits", Gaudi::DataHandle::Reader, this};
    DataHandle<Clusters>                 m_inputClusters{"inputClusters", Gaudi::DataHandle::Reader, this};
    DataHandle<TrackParametersContainer> m_outputInitialTrackParameters{"outputInitialTrackParameters",
                                                                        Gaudi::DataHandle::Writer, this};

  public:
    TrackParamVertexClusterInit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("inputVertexHits", m_inputVertexHits, "Vertex tracker hits");
      declareProperty("inputClusters", m_inputClusters, "Input clusters");
      declareProperty("outputInitialTrackParameters", m_outputInitialTrackParameters, "");
    }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
      if(!randSvc)
        return StatusCode::FAILURE;
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const Clusters*   clusters = m_inputClusters.get();
      const VertexHits* vtx_hits = m_inputVertexHits.get();
      // Create output collections
      auto init_trk_params = m_outputInitialTrackParameters.createAndPut();

      for(const auto& c : *clusters) {

        using Acts::UnitConstants::GeV;
        using Acts::UnitConstants::MeV;
        using Acts::UnitConstants::mm;
        using Acts::UnitConstants::ns;

        double p = c.energy()*GeV;
        if( p < 1.0) {
          debug() << " skipping cluster with energy " << p/GeV << " GeV" << endmsg;
          continue;
        }

        double len =  std::hypot( c.x() , c.y() , c.z() );

        // build some track cov matrix
        Acts::BoundSymMatrix cov        = Acts::BoundSymMatrix::Zero();
        cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 1.0 * mm*1.0 * mm;
        cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 1.0 * mm*1.0 * mm;
        cov(Acts::eBoundPhi, Acts::eBoundPhi)     = M_PI / 180.0;
        cov(Acts::eBoundTheta, Acts::eBoundTheta) = M_PI / 180.0;
        cov(Acts::eBoundQOverP, Acts::eBoundQOverP)     = 1.0 / (p*p);
        cov(Acts::eBoundTime, Acts::eBoundTime)         = Acts::UnitConstants::ns;

        // add all charges to the track candidate...
        init_trk_params->emplace_back(Acts::Vector4D(0 * mm, 0 * mm, 0 * mm, 0),
                                      Acts::Vector3D(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, -1,
                                      std::make_optional(cov));
        init_trk_params->emplace_back(Acts::Vector4D(0 * mm, 0 * mm, 0 * mm, 0),
                                      Acts::Vector3D(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, 1,
                                      std::make_optional(cov));
        //init_trk_params->emplace_back(Acts::Vector4D(0 * mm, 0 * mm, 0 * mm, 0),
        //                              Acts::Vector3D(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, 0,
        //                              std::make_optional(cov));
        debug() << "Invoke track finding seeded by truth particle with p = " << p/GeV  << " GeV" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(TrackParamVertexClusterInit)

} // namespace Jug::reco

