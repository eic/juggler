#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/DataHandle.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"

namespace Jug::Base {

    class MC2DummyParticle : public GaudiAlgorithm {
    public:
      DataHandle<edm4hep::MCParticleCollection> m_inputParticles{"MCParticles", Gaudi::DataHandle::Reader, this};
      DataHandle<eicd::ReconstructedParticleCollection> m_outputParticles{"DummyReconstructedParticles",
                                                                             Gaudi::DataHandle::Writer, this};
      Rndm::Numbers                                    m_gaussDist;
      Gaudi::Property<double>                          m_smearing{this, "smearing", 0.01 /* 1 percent*/};

      MC2DummyParticle(const std::string& name, ISvcLocator* svcLoc) 
        : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("inputCollection", m_inputParticles, "MCParticles");
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
        const edm4hep::MCParticleCollection* parts = m_inputParticles.get();
        // output collection
        auto& out_parts = *(m_outputParticles.createAndPut());
        for (const auto& p : *parts) {
          if (p.getGeneratorStatus() > 1) {
            if (msgLevel(MSG::DEBUG)) {
              debug() << "ignoring particle with generatorStatus = " << p.getGeneratorStatus() << endmsg;
            }
            continue;
          }

          // for now just use total momentum smearing as this is the largest effect,
          // ideally we should also smear the angles but this should be good enough
          // for now.
          const auto pvec     = p.getMomentum();
          const auto pgen     = std::hypot(pvec.x, pvec.y, pvec.z);
          const auto momentum = pgen * m_gaussDist();
          // make sure we keep energy consistent
          using MomType        = decltype(eicd::ReconstructedParticle().momentum().x);
          const MomType energy = std::sqrt(p.getEnergy() * p.getEnergy() - pgen * pgen + momentum * momentum);
          const MomType px     = p.getMomentum().x * momentum / pgen;
          const MomType py     = p.getMomentum().y * momentum / pgen;
          const MomType pz     = p.getMomentum().z * momentum / pgen;

          const MomType dpx = m_smearing.value() * px;
          const MomType dpy = m_smearing.value() * py;
          const MomType dpz = m_smearing.value() * pz;
          const MomType dE  = m_smearing.value() * energy;
          // ignore covariance for now
          // @TODO: vertex smearing
          const MomType vx = p.getVertex().x;
          const MomType vy = p.getVertex().y;
          const MomType vz = p.getVertex().z;

          auto rec_part = out_parts.create();
          rec_part.type(-1);
          rec_part.energy(energy);
          rec_part.momentum({px, py, pz});
          rec_part.referencePoint({vx, vy, vz}); // @FIXME: probably not what we want?
          rec_part.charge(p.getCharge());
          rec_part.mass(p.getMass());
          rec_part.goodnessOfPID(1); // Perfect PID
          rec_part.covMatrix({dpx, dpy, dpz, dE});
          rec_part.PDG(p.getPDG());
        }
        return StatusCode::SUCCESS;
      }
    };

    DECLARE_COMPONENT(MC2DummyParticle)

} // namespace Jug::Base

