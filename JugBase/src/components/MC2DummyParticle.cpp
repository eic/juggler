#include <algorithm>
#include <cmath>
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/GaudiTool.h"
#include "Gaudi/Algorithm.h"

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

      MC2DummyParticle(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("inputCollection", m_inputHitCollection, "mcparticles");
        declareProperty("outputCollection", m_outputHitCollection, "DummyReconstructedParticles");
      }
      StatusCode initialize() override
      {
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;
        return StatusCode::SUCCESS;
      }
      StatusCode execute() override
      {
        // input collection
        const dd4pod::Geant4ParticleCollection* parts = m_inputHitCollection.get();
        // output collection
        auto out_parts = m_outputHitCollection.createAndPut();
        for (const auto& p : *parts) {
          double momentum = std::hypot(p.psx(), p.psy(), p.psz());
          double energy   = std::hypot(momentum, p.mass());
          eic::ReconstructedParticle rec_part(p.pdgID(), energy, {p.psx(),p.psy(),p.psz()}, (double)p.charge(), p.mass());
          out_parts->push_back(rec_part);
        }
        return StatusCode::SUCCESS;
      }
    };

    DECLARE_COMPONENT(MC2DummyParticle)

  } // namespace Base
} // namespace Gaudi

