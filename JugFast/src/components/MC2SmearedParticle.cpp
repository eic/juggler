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
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"

namespace Jug::Fast {

class MC2SmearedParticle : public GaudiAlgorithm, AlgorithmIDMixin<int32_t> {
public:
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"MCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"SmearedReconstructedParticles",
                                                                         Gaudi::DataHandle::Writer, this};
  Rndm::Numbers m_gaussDist;
  Gaudi::Property<double> m_smearing{this, "smearing", 0.01 /* 1 percent*/};

  const int32_t kMonteCarloSource{uniqueID<int32_t>("MCParticles")};

  MC2SmearedParticle(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputParticles", m_inputMCParticles, "MCParticles");
    declareProperty("outputParticles", m_outputParticles, "SmearedReconstructedParticles");
  }
  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;

    IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
    StatusCode sc        = m_gaussDist.initialize(randSvc, Rndm::Gauss(1.0, m_smearing.value()));
    if (!sc.isSuccess()) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
  StatusCode execute() override {
    // input collection
    auto parts = m_inputMCParticles.get();
    // output collection
    auto& out_parts = *(m_outputParticles.createAndPut());
    int ID         = 0;
    for (const auto& p : *parts) {
      if (p.getGeneratorStatus() > 1) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "ignoring particle with getGeneratorStatus = " << p.getGeneratorStatus() << endmsg;
        }
        continue;
      }

      // for now just use total momentum smearing as this is the largest effect,
      // ideally we should also smear the angles but this should be good enough
      // for now.
      const auto pvec      = p.getMomentum();
      const double pgen    = std::hypot(pvec.x, pvec.y, pvec.z);
      const float momentum = pgen * m_gaussDist();
      const float energy   = p.getEnergy();
      const float px       = p.getMomentum().x * momentum / pgen;
      const float py       = p.getMomentum().y * momentum / pgen;
      const float pz       = p.getMomentum().z * momentum / pgen;
      // @TODO: vertex smearing
      const float vx = p.getVertex().x;
      const float vy = p.getVertex().y;
      const float vz = p.getVertex().z;

      eic::VectorXYZ psmear{px, py, pz};

      auto rec_part = out_parts.create();
      rec_part.ID({ID++, algorithmID()});
      rec_part.p({px, py, pz});
      rec_part.v({vx, vy, vz});
      rec_part.time(p.getTime());
      rec_part.pid(p.getPDG());
      rec_part.status(p.getGeneratorStatus());
      rec_part.charge(p.getCharge());
      rec_part.weight(1.);
      rec_part.direction({psmear.theta(), psmear.phi()});
      rec_part.momentum(momentum);
      rec_part.energy(energy);
      rec_part.mass(p.getMass());
      // FIXME disabled mcID link
      //rec_part.mcID({p.id(), kMonteCarloSource});
    }
    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(MC2SmearedParticle)

} // namespace Jug::Fast
