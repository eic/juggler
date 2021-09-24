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

//#include "Acts/Definitions/Units.hpp"
//#include "Acts/Definitions/Common.hpp"
//#include "Acts/EventData/Charge.hpp"

#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "Math/Vector3D.h"


namespace Jug::Reco {

  /** Track seeding using MC truth.
   *
   *  \note "Seeding" algorithms are required to output a eic::TrackParametersCollection, as opposed to the legacy "init"
   *  algorithms, such as  TrackParamTruthInit.
   *
   *  \ingroup tracking
   */
  class TruthTrackSeeding : public GaudiAlgorithm {
  public:
    DataHandle<dd4pod::Geant4ParticleCollection> m_inputMCParticles{"inputMCParticles", Gaudi::DataHandle::Reader,
                                                                    this};
    DataHandle<eic::TrackParametersCollection> m_outputTrackParameters{"outputTrackParameters",
                                                                       Gaudi::DataHandle::Writer, this};
    SmartIF<IParticleSvc> m_pidSvc;

  public:
    TruthTrackSeeding(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
      declareProperty("inputMCParticles", m_inputMCParticles, "mcparticle truth data from npsim");
      declareProperty("outputTrackParameters", m_outputTrackParameters, "Output initial track parameters");
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
      auto init_trk_params = m_outputTrackParameters.createAndPut();
      int ID               = 0;

      for(const auto& part : *mcparts) {

        // genStatus = 1 means thrown G4Primary 
        if(part.genStatus() != 1 ) {
          continue;
        }

        const double p = part.ps().mag();

        // get the particle charge
        // note that we cannot trust the mcparticles charge, as DD4hep
        // sets this value to zero! let's lookup by PDGID instead
        const double charge = m_pidSvc->particle(part.pdgID()).charge;
        if (abs(charge) < std::numeric_limits<double>::epsilon()) {
          continue;
        }

        float q_over_p = charge / p;

        // TrackParameters params(ID, loc, eic::FloatPair locError, eic::Direction direction, eic::Direction
        // directionError, float qOverP, float qOverPError, float time, float timeError);
        eic::TrackParameters params(eic::Index(ID++), {float(0.0), float(0.0)}, {float(0.1), float(0.1)},
                                    {float(part.ps().phi()), float(part.ps().theta())}, {float(0.1), float(0.1)},
                                    q_over_p, float(q_over_p * 0.05), float(part.time()), float(1.0), float(charge));

        ////// Construct a perigee surface as the target surface
        //auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        //    Acts::Vector3{part.vs().x * mm, part.vs().y * mm, part.vs().z * mm});

        init_trk_params->push_back(params);

        if (msgLevel(MSG::DEBUG)) {
          debug() << "Invoke track finding seeded by truth particle with p = " << p  << " GeV" << endmsg;
          debug() << "                                              charge = " << charge << endmsg;
          debug() << "                                                 q/p = " << charge / p << endmsg;
        }
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(TruthTrackSeeding)

} // namespace Jug::reco

