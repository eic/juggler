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
#include "dd4pod/Geant4ParticleCollection.h"
#include "Math/Vector3D.h"
#include "Acts/Surfaces/PerigeeSurface.hpp"

#include "TDatabasePDG.h"

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
   */
  class TrackParamTruthInit : public GaudiAlgorithm {
  public:
    DataHandle<dd4pod::Geant4ParticleCollection> m_inputMCParticles{"inputMCParticles", Gaudi::DataHandle::Reader,
                                                                    this};
    DataHandle<TrackParametersContainer>         m_outputInitialTrackParameters{"outputInitialTrackParameters",
                                                                        Gaudi::DataHandle::Writer, this};
    TDatabasePDG pdgData;

  public:
    TrackParamTruthInit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("inputMCParticles", m_inputMCParticles, "");
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
      const dd4pod::Geant4ParticleCollection* mcparts = m_inputMCParticles.get();
      // Create output collections
      auto init_trk_params = m_outputInitialTrackParameters.createAndPut();

      for(const auto& part : *mcparts) {

        // genStatus = 1 means thrown G4Primary
        if(part.genStatus() != 1 ) {
          continue;
        }
        using Acts::UnitConstants::GeV;
        using Acts::UnitConstants::MeV;
        using Acts::UnitConstants::mm;
        using Acts::UnitConstants::um;
        using Acts::UnitConstants::ns;

        double p = std::hypot( part.psx() * GeV, part.psy() * GeV, part.psz() * GeV);
        ROOT::Math::XYZVector  momentum(part.psx() * GeV, part.psy() * GeV, part.psz() * GeV);

        // build some track cov matrix
        Acts::BoundSymMatrix cov        = Acts::BoundSymMatrix::Zero();
        cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 1000*um*1000*um;
        cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 1000*um*1000*um;
        cov(Acts::eBoundPhi, Acts::eBoundPhi)     = 0.05*0.05;
        cov(Acts::eBoundTheta, Acts::eBoundTheta) = 0.01*0.01;
        cov(Acts::eBoundQOverP, Acts::eBoundQOverP)     = (0.1*0.1) / (GeV*GeV);
        cov(Acts::eBoundTime, Acts::eBoundTime)         = 10.0e9*ns*10.0e9*ns;

        Acts::BoundVector  params;
        params(Acts::eBoundLoc0)   = 0.0 * mm ;  // cylinder radius
        params(Acts::eBoundLoc1)   = 0.0 * mm ; // cylinder length
        params(Acts::eBoundPhi)    = momentum.Phi();
        params(Acts::eBoundTheta)  = momentum.Theta();
        params(Acts::eBoundQOverP) = -1/p;
        params(Acts::eBoundTime)   = part.time() * ns;
        // get the particle charge
        double charge = pdgData.GetParticle(part.pdgID()) -> Charge();

        //// Construct a perigee surface as the target surface
        auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
            Acts::Vector3{part.vsx() * mm, part.vsy() * mm, part.vsz() * mm});

        //params(Acts::eBoundQOverP) = charge/p;
        init_trk_params->push_back({pSurface, params, charge, cov});
        // std::make_optional(std::move(cov))

        debug() << "Invoke track finding seeded by truth particle with p = " << p/GeV  << " GeV" << endmsg;
        //Acts::BoundMatrix cov           = Acts::BoundMatrix::Zero();
        //cov(Acts::eLOC_0, Acts::eLOC_0) = ahit.covMatrix(0)*ahit.covMatrix(0);
        //cov(Acts::eLOC_1, Acts::eLOC_1) = ahit.covMatrix(1)*ahit.covMatrix(1);
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(TrackParamTruthInit)

} // namespace Jug::reco

