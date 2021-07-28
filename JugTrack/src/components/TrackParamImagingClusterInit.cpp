#include <cmath>
// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugTrack/Track.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"

#include "eicd/TrackerHitCollection.h"
#include "eicd/ImagingClusterCollection.h"

#include "Math/Vector3D.h"
#include "Acts/Surfaces/PerigeeSurface.hpp"

  ///// (Reconstructed) track parameters e.g. close to the vertex.
  //using TrackParameters = Acts::CurvilinearTrackParameters;

  ///// Container of reconstructed track states for multiple tracks.
  //using TrackParametersContainer = std::vector<TrackParameters>;

  ///// MultiTrajectory definition
  //using Trajectory = Acts::MultiTrajectory<SourceLink>;

  ///// Container for the truth fitting/finding track(s)
  //using TrajectoryContainer = std::vector<SimMultiTrajectory>;

namespace Jug::Reco {

  /** Initial Track parameters from MC truth.
   *
   *  TrackParmetersContainer
   *
   *  \ingroup track
   */
  class TrackParamImagingClusterInit : public GaudiAlgorithm {
  public:
    using ImagingClusters =  eic::ImagingClusterCollection;

    DataHandle<ImagingClusters>          m_inputClusters{"inputClusters", Gaudi::DataHandle::Reader, this};
    DataHandle<TrackParametersContainer> m_outputInitialTrackParameters{"outputInitialTrackParameters",
                                                                        Gaudi::DataHandle::Writer, this};

  public:
    TrackParamImagingClusterInit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
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
      const eic::ImagingClusterCollection* clusters = m_inputClusters.get();
      // Create output collections
      auto init_trk_params = m_outputInitialTrackParameters.createAndPut();

      for(const auto& c : *clusters) {

        using Acts::UnitConstants::GeV;
        using Acts::UnitConstants::MeV;
        using Acts::UnitConstants::mm;
        using Acts::UnitConstants::ns;

        double p = c.energy()*GeV;
        if( p < 0.1*GeV) {
          continue;
        }
        double len =  std::hypot( c.x() , c.y() , c.z() );
        ROOT::Math::XYZVector  momentum(c.x() * p / len, c.y() * p / len, c.z() * p / len);

        // build some track cov matrix
        //Acts::BoundSymMatrix cov        = Acts::BoundSymMatrix::Zero();
        //cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 1.0 * mm*1.0 * mm;
        //cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 1.0 * mm*1.0 * mm;
        //cov(Acts::eBoundPhi, Acts::eBoundPhi)     = M_PI / 180.0;
        //cov(Acts::eBoundTheta, Acts::eBoundTheta) = M_PI / 180.0;
        //cov(Acts::eBoundQOverP, Acts::eBoundQOverP)     = 1.0 / (p*p);
        //cov(Acts::eBoundTime, Acts::eBoundTime)         = Acts::UnitConstants::ns;

        Acts::BoundVector  params;
        params(Acts::eBoundLoc0)   = 0.0 * mm ;
        params(Acts::eBoundLoc1)   = 0.0 * mm ;
        params(Acts::eBoundPhi)    = momentum.Phi();
        params(Acts::eBoundTheta)  = momentum.Theta();
        params(Acts::eBoundQOverP) = 1/p;
        params(Acts::eBoundTime)   = 0 * ns;

        auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0,0,0});

        debug() << "Invoke track finding seeded by truth particle with p = " << p/GeV  << " GeV" << endmsg;

        // add both charges to the track candidate...
        init_trk_params->push_back({pSurface, params,  1});

        Acts::BoundVector  params2;
        params2(Acts::eBoundLoc0)   = 0.0 * mm ;
        params2(Acts::eBoundLoc1)   = 0.0 * mm ;
        params2(Acts::eBoundPhi)    = momentum.Phi();
        params2(Acts::eBoundTheta)  = momentum.Theta();
        params2(Acts::eBoundQOverP) = -1/p;
        params2(Acts::eBoundTime)   = 0 * ns;
        init_trk_params->push_back({pSurface, params2, -1});

        // acts v1.2.0:
        //init_trk_params->emplace_back(Acts::Vector4(0 * mm, 0 * mm, 0 * mm, 0),
        //                              Acts::Vector3(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, -1,
        //                              std::make_optional(cov));
        //debug() << init_trk_params->back() << endmsg;
        //init_trk_params->emplace_back(Acts::Vector4(0 * mm, 0 * mm, 0 * mm, 0),
        //                              Acts::Vector3(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, 1,
        //                              std::make_optional(cov));
        ////debug() << init_trk_params->back() << endmsg;
        //init_trk_params->emplace_back(Acts::Vector4(0 * mm, 0 * mm, 0 * mm, 0),
        //                              Acts::Vector3(c.x() * p / len, c.y() * p / len, c.z() * p / len), p, 0,
        //                              std::make_optional(cov));
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(TrackParamImagingClusterInit)

} // namespace Jug::reco

