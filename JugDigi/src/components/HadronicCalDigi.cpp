#include <algorithm>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/GaudiTool.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitData.h"
#include "dd4pod/CalorimeterHitCollection.h"

namespace Jug {
  namespace Digi {

   class HadronicCalDigi : public GaudiAlgorithm {
   public:
     Gaudi::Property<double>                      m_energyResolution{this, "energyResolution", 50/*percent/sqrt(E)*/};
     Rndm::Numbers                                m_gaussDist;
     DataHandle<dd4pod::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
     DataHandle<eic::RawCalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
   public:

    HadronicCalDigi(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputHitCollection", m_inputHitCollection,"");
          declareProperty("outputHitCollection", m_outputHitCollection, "");
        }

    StatusCode initialize() override {
      IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
      StatusCode   sc      = m_gaussDist.initialize( randSvc, Rndm::Gauss(1.0, m_energyResolution.value()/100.0));
      if (!sc.isSuccess()) {
        return StatusCode::FAILURE;
      }
      if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;
      //f_counter = m_starting_value.value();
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const dd4pod::CalorimeterHitCollection* simhits = m_inputHitCollection.get();
      // Create output collections
      auto rawhits = m_outputHitCollection.createAndPut();
      eic::RawCalorimeterHitCollection* rawHitCollection = new eic::RawCalorimeterHitCollection();
      for(const auto& ahit : *simhits) {
        //std::cout << ahit << "\n";
        // here 1000 is arbitrary
        eic::RawCalorimeterHit rawhit((long long)ahit.cellID(), (long long)ahit.cellID(),
                        (long long)ahit.energyDeposit()*1000, 0);
        rawhits->push_back(rawhit);
      }
      return StatusCode::SUCCESS;
    }

  };
  DECLARE_COMPONENT(HadronicCalDigi)

  } // namespace Examples
} // namespace Gaudi

