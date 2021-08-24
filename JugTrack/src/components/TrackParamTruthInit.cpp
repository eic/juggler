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
#include "JugBase/IParticleSvc.h"
#include "JugTrack/Track.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"

#include "eicd/TrackerHitCollection.h"
#include "dd4pod/Geant4ParticleCollection.h"
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
   *  \ingroup tracking
   */
  class TrackParamTruthInit : public GaudiAlgorithm {
  public:
    DataHandle<dd4pod::Geant4ParticleCollection> m_inputMCParticles{"inputMCParticles", Gaudi::DataHandle::Reader,
                                                                    this};
    DataHandle<TrackParametersContainer>         m_outputInitialTrackParameters{"outputInitialTrackParameters",
                                                                        Gaudi::DataHandle::Writer, this};

    SmartIF<IParticleSvc> m_pidSvc;

  public:
    TrackParamTruthInit(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("inputMCParticles", m_inputMCParticles, "");
      declareProperty("outputInitialTrackParameters", m_outputInitialTrackParameters, "");
    }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
      if (!randSvc) {
        return StatusCode::FAILURE;
      }
      m_pidSvc = service("ParticleSvc");
      if (!m_pidSvc) {
        error() << "Unable to locate Particle Service. "
                << "Make sure you have ParticleSvc in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }
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

        const double p = part.ps().mag() * GeV;

        // get the particle charge
        // note that we cannot trust the mcparticles charge, as DD4hep
        // sets this value to zero! let's lookup by PDGID instead
        const double charge = m_pidSvc->particle(part.pdgID()).charge;

        // build some track cov matrix
        Acts::BoundSymMatrix cov                    = Acts::BoundSymMatrix::Zero();
        cov(Acts::eBoundLoc0, Acts::eBoundLoc0)     = 1000*um*1000*um;
        cov(Acts::eBoundLoc1, Acts::eBoundLoc1)     = 1000*um*1000*um;
        cov(Acts::eBoundPhi, Acts::eBoundPhi)       = 0.05*0.05;
        cov(Acts::eBoundTheta, Acts::eBoundTheta)   = 0.01*0.01;
        cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = (0.1*0.1) / (GeV*GeV);
        cov(Acts::eBoundTime, Acts::eBoundTime)     = 10.0e9*ns*10.0e9*ns;

        Acts::BoundVector  params;
        params(Acts::eBoundLoc0)   = 0.0 * mm ;  // cylinder radius
        params(Acts::eBoundLoc1)   = 0.0 * mm ; // cylinder length
        params(Acts::eBoundPhi)    = part.ps().phi();
        params(Acts::eBoundTheta)  = part.ps().theta();
        params(Acts::eBoundQOverP) = charge / p;
        params(Acts::eBoundTime)   = part.time() * ns;

        //// Construct a perigee surface as the target surface
        auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
            Acts::Vector3{part.vs().x * mm, part.vs().y * mm, part.vs().z * mm});

        //params(Acts::eBoundQOverP) = charge/p;
        init_trk_params->push_back({pSurface, params, charge,cov});
        // std::make_optional(std::move(cov))

        if (msgLevel(MSG::DEBUG)) {
          debug() << "Invoke track finding seeded by truth particle with p = " << p / GeV << " GeV" << endmsg;
          debug() << "                                              charge = " << charge << endmsg;
          debug() << "                                                 q/p = " << charge / p << endmsg;
        }
        //Acts::BoundMatrix cov           = Acts::BoundMatrix::Zero();
        //cov(Acts::eLOC_0, Acts::eLOC_0) = ahit.covMatrix(0)*ahit.covMatrix(0);
        //cov(Acts::eLOC_1, Acts::eLOC_1) = ahit.covMatrix(1)*ahit.covMatrix(1);
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(TrackParamTruthInit)

} // namespace Jug::reco

