#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/DataHandle.h"
#include "JugBase/Utilities/UniqueID.hpp"

// Event Model related classes
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"

namespace Jug {
  namespace Base {

    class MC2DummyParticle : public GaudiAlgorithm {
  private:
    // Unique identifier for this hit type, based on the algorithm name
    using PartClassificationType = decltype(eic::ReconstructedParticleData().type);
    const PartClassificationType m_type;
    public:
      DataHandle<dd4pod::Geant4ParticleCollection> m_inputHitCollection{"mcparticles", Gaudi::DataHandle::Reader, this};
      DataHandle<eic::ReconstructedParticleCollection> m_outputHitCollection{"DummyReconstructedParticles",
                                                                             Gaudi::DataHandle::Writer, this};
      Rndm::Numbers                                    m_gaussDist;
      Gaudi::Property<double>                          m_smearing{this, "smearing", 0.01 /* 1 percent*/};

      MC2DummyParticle(const std::string& name, ISvcLocator* svcLoc) 
        : GaudiAlgorithm(name, svcLoc)
        , m_type{uniqueID<PartClassificationType>(name)}
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
        int ID = 0;
        for (const auto& p : *parts) {
          if (p.genStatus() != 1) {
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

          eic::ReconstructedParticle rec_part{
            ID++,                             // Unique index
            {px, py, pz},                     // 3-momentum [GeV]
            {vx, vy, vz},                     // Vertex [mm]
            static_cast<float>(p.time()),     // time [ns]
            p.pdgID(),                        // PDG type
            static_cast<int16_t>(p.status()), // Status
            static_cast<int16_t>(p.charge()), // Charge
            m_type,                           // Algorithm type
            momentum,                         // 3-momentum magnitude [GeV]
            energy,                           // energy [GeV]
            static_cast<float>(p.mass()),     // mass [GeV]
            1.};                              // particle weight

          out_parts->push_back(rec_part);
        }
        return StatusCode::SUCCESS;
      }
    };

    DECLARE_COMPONENT(MC2DummyParticle)

  } // namespace Base
} // namespace Jug

