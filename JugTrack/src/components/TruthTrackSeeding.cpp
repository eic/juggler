// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten

#include <cmath>
// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>
#include "JugBase/IParticleSvc.h"
#include "ActsExamples/EventData/Track.hpp"

//#include "Acts/Definitions/Units.hpp"
//#include "Acts/Definitions/Common.hpp"
//#include "Acts/EventData/Charge.hpp"

#include "edm4eic/EDM4eicVersion.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/TrackParametersCollection.h"
#include "Math/Vector3D.h"


namespace Jug::Reco {

  /** Track seeding using MC truth.
   *
   *  \note "Seeding" algorithms are required to output a edm4eic::TrackParametersCollection, as opposed to the legacy "init"
   *  algorithms, such as  TrackParamTruthInit.
   *
   *  \ingroup tracking
   */
  class TruthTrackSeeding : public Gaudi::Algorithm {
  private:
    mutable DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"inputMCParticles", Gaudi::DataHandle::Reader,
                                                                    this};
    mutable DataHandle<edm4eic::TrackParametersCollection> m_outputTrackParameters{"outputTrackParameters",
                                                                       Gaudi::DataHandle::Writer, this};
    SmartIF<IParticleSvc> m_pidSvc;

  public:
    TruthTrackSeeding(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
      declareProperty("inputMCParticles", m_inputMCParticles, "mcparticle truth data from npsim");
      declareProperty("outputTrackParameters", m_outputTrackParameters, "Output initial track parameters");
    }

    StatusCode initialize() override {
      if (Gaudi::Algorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      IRndmGenSvc* randSvc = Gaudi::svcLocator()->service<IRndmGenSvc>("RndmGenSvc", true);
      if (randSvc == nullptr) {
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

    StatusCode execute(const EventContext&) const override {
      // input collection
      const auto* const mcparts = m_inputMCParticles.get();
      // Create output collections
      auto* init_trk_params = m_outputTrackParameters.createAndPut();

      for(const auto& part : *mcparts) {

        // getGeneratorStatus = 1 means thrown G4Primary 
        if(part.getGeneratorStatus() != 1 ) {
          continue;
        }

        const auto& pvec = part.getMomentum();
        const auto p = std::hypot(pvec.x, pvec.y, pvec.z);
        const auto phi = std::atan2(pvec.x, pvec.y);
        const auto theta = std::atan2(std::hypot(pvec.x, pvec.y), pvec.z);

        // get the particle charge
        // note that we cannot trust the mcparticles charge, as DD4hep
        // sets this value to zero! let's lookup by PDGID instead
        const auto charge = static_cast<float>(m_pidSvc->particle(part.getPDG()).charge);
        if (abs(charge) < std::numeric_limits<double>::epsilon()) {
          continue;
        }

        const auto q_over_p = charge / p;

        auto params = init_trk_params->create();
        params.setType(-1);                  // type --> seed (-1)
        params.setLoc({0.0F, 0.0F});         // location on surface
        params.setTheta(theta);              // theta (rad)
        params.setPhi(phi);                  // phi  (rad)
        params.setQOverP(q_over_p * .05F);   // Q/P (e/GeV)
        params.setTime(part.getTime());      // Time (ns)
        #if EDM4EIC_VERSION_MAJOR >= 5
          edm4eic::Cov6f cov;
          cov(0,0) = 0.1; // loc0
          cov(0,1) = 0.1; // loc0/loc1
          cov(1,1) = 0.1; // loc1
          cov(2,2) = 0.1; // theta
          cov(3,3) = 0.1; // phi
          cov(4,4) = 0.1; // qOverP
          cov(5,5) = 0.1; // time
          params.setCovariance(cov);
        #else
          params.setLocError({0.1, 0.1, 0.1}); // Covariance on location
          params.setMomentumError({0.1, 0.1, 0.1}); // Covariance on theta/phi/Q/P
          params.setTimeError(0.1);            // Error on time
          params.setCharge(charge);            // Charge
        #endif

        ////// Construct a perigee surface as the target surface
        //auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        //    Acts::Vector3{part.getVertex().x * mm, part.getVertex().y * mm, part.getVertex().z * mm});

        if (msgLevel(MSG::DEBUG)) {
          debug() << "Invoke track finding seeded by truth particle with p = " << p  << " GeV" << endmsg;
          debug() << "                                              charge = " << charge << endmsg;
          debug() << "                                                 q/p = " << charge / p << endmsg;
        }
      }
      return StatusCode::SUCCESS;
    }
  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(TruthTrackSeeding)

} // namespace Jug::reco
