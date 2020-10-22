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
#include "dd4pod/Geant4ParticleCollection.h"


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
        using Acts::UnitConstants::MeV;
        using Acts::UnitConstants::GeV;
        using Acts::UnitConstants::mm;

        double p = std::hypot( part.psx() * MeV, part.psy() * MeV, part.psz() * MeV);

        // build some track cov matrix
        Acts::BoundSymMatrix cov        = Acts::BoundSymMatrix::Zero();
        cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 1.0 * mm*1.0 * mm;
        cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 1.0 * mm*1.0 * mm;
        cov(Acts::eBoundPhi, Acts::eBoundPhi)     = M_PI / 180.0;
        cov(Acts::eBoundTheta, Acts::eBoundTheta) = M_PI / 180.0;
        cov(Acts::eBoundQOverP, Acts::eBoundQOverP)     = 0.98 / (p*p);
        cov(Acts::eBoundTime, Acts::eBoundTime)         = Acts::UnitConstants::ns;

        init_trk_params->emplace_back(Acts::Vector4D(part.vsx() * mm, part.vsy() * mm, part.vsz() * mm, part.time() * Acts::UnitConstants::ns),
                                      Acts::Vector3D(part.psx() * MeV, part.psy() * MeV, part.psz() * MeV),
                                      p+200*MeV,
                                      ((part.pdgID() > 0) ? -1 : 1),
                                      std::make_optional(std::move(cov))
                                      );
        //part .charge()

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

