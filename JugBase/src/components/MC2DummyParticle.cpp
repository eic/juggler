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
      //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
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
          out_parts->push_back({p.pdgID(), energy, {p.psx(),p.psy(),p.psz()}, p.charge(), p.mass()});
        }
        return StatusCode::SUCCESS;
      }

      DataHandle<dd4pod::Geant4ParticleCollection>     m_inputHitCollection{"mcparticles", Gaudi::DataHandle::Reader, this};
      DataHandle<eic::ReconstructedParticleCollection> m_outputHitCollection{"DummyReconstructedParticles", Gaudi::DataHandle::Writer, this};
    };

    DECLARE_COMPONENT(MC2DummyParticle)

  } // namespace Base
} // namespace Gaudi

