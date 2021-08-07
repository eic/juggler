#include <algorithm>
#include <cmath>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitData.h"
#include "dd4pod/CalorimeterHitCollection.h"

namespace Jug {
  namespace Digi {

    /** Electromagnetic calorimeter digitization.
     *
     * \ingroup digi
     */
    class EMCalorimeterDigi : public GaudiAlgorithm {
    public:
      using SimHit = dd4pod::CalorimeterHitCollection;
      using RawHit = eic::RawCalorimeterHitCollection;
      DataHandle<SimHit> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
      DataHandle<RawHit> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};

      Gaudi::Property<double> m_energyResolution{this, "energyResolution", 0.05 /* 5 percent*/};
      Rndm::Numbers           m_gaussDist;

      EMCalorimeterDigi(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("inputHitCollection", m_inputHitCollection, "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
      }
      StatusCode initialize() override
      {
        warning() << "Deprecated algorithm for digi/reco, use Jug::Digi::CalorimeterHitDigi"
                     "and Jug::Reco::CalorimeterHitReco instead" << endmsg;
        IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        StatusCode   sc      = m_gaussDist.initialize(randSvc, Rndm::Gauss(0.0, m_energyResolution.value()));
        if (!sc.isSuccess()) {
          return StatusCode::FAILURE;
        }
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;
        return StatusCode::SUCCESS;
      }
      StatusCode execute() override
      {
        // input collection
        const dd4pod::CalorimeterHitCollection* simhits = m_inputHitCollection.get();
        // Create output collections
        auto                              rawhits          = m_outputHitCollection.createAndPut();
        eic::RawCalorimeterHitCollection* rawHitCollection = new eic::RawCalorimeterHitCollection();
        int nhits = 0;
        for (const auto& ahit : *simhits) {
          // std::cout << ahit << "\n";
          double sqrtE = std::sqrt(ahit.energyDeposit()) ;
          double aterm = m_gaussDist()*sqrtE;
          eic::RawCalorimeterHit rawhit((long long)ahit.cellID(),
                                        std::llround((ahit.energyDeposit() + aterm) * 1e6),
                                        ahit.truth().time * 1e6, nhits++);
          rawhits->push_back(rawhit);
        }
        return StatusCode::SUCCESS;
      }

    };
    DECLARE_COMPONENT(EMCalorimeterDigi)
  } // namespace Digi
} // namespace Jug
