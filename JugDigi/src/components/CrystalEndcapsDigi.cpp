#include <algorithm>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Property.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitData.h"
#include "dd4pod/CalorimeterHitCollection.h"

namespace Jug {
  namespace Digi {
  
    /** Crystal Endcaps Calorimeter detector digitization.
     *
     *
     */
   class CrystalEndcapsDigi : public GaudiAlgorithm {
   public:
  
    Gaudi::Property<double>      m_energyResolution{this, "energyResolution", 0.02};  // 2%sqrt(E)
    Rndm::Numbers m_gaussDist;
    DataHandle<dd4pod::CalorimeterHitCollection> m_inputHitCollection{ "inputHitCollection",  Gaudi::DataHandle::Reader, this};
    DataHandle<eic::RawCalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};

    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    CrystalEndcapsDigi(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputHitCollection", m_inputHitCollection,"");
          declareProperty("outputHitCollection", m_outputHitCollection, "");
        }
    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
      StatusCode   sc      = m_gaussDist.initialize( randSvc, Rndm::Gauss(0.0, m_energyResolution.value()));
      if (!sc.isSuccess()) {
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }
    StatusCode execute() override {
      // input collections
      const dd4pod::CalorimeterHitCollection* simhits = m_inputHitCollection.get();
      // Create output collections
      auto rawhits = m_outputHitCollection.createAndPut();
      eic::RawCalorimeterHitCollection* rawHitCollection = new eic::RawCalorimeterHitCollection();
      for(const auto& ahit : *simhits){
	   	eic::RawCalorimeterHit rawhit((long long)ahit.cellID(), (long long)ahit.cellID(), 
			(long long)(ahit.energyDeposit() + m_gaussDist*sqrt(ahit.energyDeposit())) * 100.0, (double)ahit.truth().time);
          	rawhits->push_back(rawhit);
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(CrystalEndcapsDigi)
  } // namespace Digi
} // namespace Jug

