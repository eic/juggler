#include <algorithm>
#include <cmath>
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/GaudiTool.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/RndmGenerators.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"


namespace Jug {
  namespace Base {

    class MC2DummyParticle : public GaudiAlgorithm {
    public:
      DataHandle<dd4pod::Geant4ParticleCollection>     m_inputHitCollection{"mcparticles", Gaudi::DataHandle::Reader, this};
      DataHandle<eic::ReconstructedParticleCollection> m_outputHitCollection{"DummyReconstructedParticles", Gaudi::DataHandle::Writer, this};
      Rndm::Numbers           m_gaussDist;
      Gaudi::Property<double> m_smearing{this, "smearing", 0.01 /* 1 percent*/};

      MC2DummyParticle(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("inputCollection", m_inputHitCollection, "mcparticles");
        declareProperty("outputCollection", m_outputHitCollection, "DummyReconstructedParticles");
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
        const dd4pod::Geant4ParticleCollection* parts = m_inputHitCollection.get();
        // output collection
        auto out_parts = m_outputHitCollection.createAndPut();
        for (const auto& p : *parts) {
          if (p.genStatus() != 1) {
            continue;
          }
          
          double momentum = std::hypot(p.psx(), p.psy(), p.psz());
          double energy   = std::hypot(momentum, p.mass());
          eic::ReconstructedParticle rec_part(p.pdgID(), energy*m_gaussDist(), {p.psx()*m_gaussDist(),p.psy()*m_gaussDist(),p.psz()*m_gaussDist()}, (double)p.charge(), p.mass());
          out_parts->push_back(rec_part);
        }
        return StatusCode::SUCCESS;
      }
    };

    DECLARE_COMPONENT(MC2DummyParticle)

  } // namespace Base
} // namespace Gaudi

