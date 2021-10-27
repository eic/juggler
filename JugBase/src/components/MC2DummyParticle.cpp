#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"

namespace Jug::Base {

    class MC2DummyParticle : public GaudiAlgorithm, AlgorithmIDMixin<int32_t> {
    public:
      DataHandle<dd4pod::Geant4ParticleCollection> m_inputParticles{"mcparticles", Gaudi::DataHandle::Reader, this};
      DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"DummyReconstructedParticles",
                                                                             Gaudi::DataHandle::Writer, this};
      Rndm::Numbers                                    m_gaussDist;
      Gaudi::Property<double>                          m_smearing{this, "smearing", 0.01 /* 1 percent*/};

      const int32_t kMonteCarloSource{uniqueID<int32_t>("mcparticles")};

      MC2DummyParticle(const std::string& name, ISvcLocator* svcLoc) 
        : GaudiAlgorithm(name, svcLoc)
        , AlgorithmIDMixin(name, info())
      {
        declareProperty("inputCollection", m_inputParticles, "mcparticles");
        declareProperty("outputCollection", m_outputParticles, "DummyReconstructedParticles");
      }
      StatusCode initialize() override
      {
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;

        IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        StatusCode   sc      = m_gaussDist.initialize(randSvc, Rndm::Gauss(1.0, m_smearing.value()));
        if (!sc.isSuccess()) {
          return StatusCode::FAILURE;
        }
        return StatusCode::SUCCESS;
      }
      StatusCode execute() override
      {
        // input collection
        const dd4pod::Geant4ParticleCollection* parts = m_inputParticles.get();
        // output collection
        auto& out_parts = *(m_outputParticles.createAndPut());
        int ID = 0;
        for (const auto& p : *parts) {
          if (p.genStatus() > 1) {
            if (msgLevel(MSG::DEBUG)) {
              debug() << "ignoring particle with genStatus = " << p.genStatus() << endmsg;
            }
            continue;
          }

          // for now just use total momentum smearing as this is the largest effect,
          // ideally we should also smear the angles but this should be good enough
          // for now.
          const double pgen     = p.ps().mag();
          const float momentum = pgen * m_gaussDist();
          const float energy   = p.energy();
          const float px       = p.ps().x * momentum / pgen;
          const float py       = p.ps().y * momentum / pgen;
          const float pz       = p.ps().z * momentum / pgen;
          // @TODO: vertex smearing
          const float vx       = p.vs().x;
          const float vy       = p.vs().y;
          const float vz       = p.vs().z;

          eic::VectorXYZ psmear{px, py, pz};

          auto rec_part = out_parts.create();
          rec_part.ID({ID++, algorithmID()});
          rec_part.p(psmear);
          rec_part.v({vx, vy, vz});
          rec_part.time(static_cast<float>(p.time()));
          rec_part.pid(p.pdgID());
          rec_part.status(p.genStatus());
          rec_part.charge(p.charge());
          rec_part.weight(1.);
          rec_part.direction({psmear.theta(), psmear.phi()});
          rec_part.momentum(momentum);
          rec_part.energy(energy);
          rec_part.mass(p.mass());
          rec_part.mcID({p.ID(), kMonteCarloSource});
        }
        return StatusCode::SUCCESS;
      }
    };

    DECLARE_COMPONENT(MC2DummyParticle)

} // namespace Jug::Base

