#include <algorithm>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
//#include "GaudiExamples/MyTrack.h"
#include "eicd/RawCalorimeterHitCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"

namespace Jug {
  namespace Digi {

  //struct ExampleCaloDigi final
  //    : Gaudi::Functional::Transformer<eic::RawCalorimeterHitCollection(
  //          const dd4pod::CalorimeterHitCollection&)> {

  //  ExampleCaloDigi(const std::string& name, ISvcLocator* pSvc)
  //      : Transformer(name, pSvc, {KeyValue("InputData", {"FAEC_ShHits"})},
  //                    KeyValue("OutputData", {"ForwardPreshowerHits"})) {}

  //  eic::RawCalorimeterHitCollection
  //  operator()(const dd4pod::CalorimeterHitCollection& in_hits) const override {
  //    eic::RawCalorimeterHitCollection out_hits;
  //    for (auto i = in_hits.begin(), end = in_hits.end(); i != end; ++i) {
  //      out_hits.create((long long)i->cellID(), (long long)i->cellID(),
  //                      (long long)i->energyDeposit() * 100, 0);
  //    }
  //    return out_hits;
  //  }
  //};
  //DECLARE_COMPONENT(ExampleCaloDigi)
  
   class ExampleCaloDigi : public GaudiAlgorithm {
   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ExampleCaloDigi(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputHitCollection", m_inputHitCollection,"");
          declareProperty("outputHitCollection", m_outputHitCollection, "");
        }
    StatusCode initialize() override {
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
        eic::RawCalorimeterHit rawhit((long long)ahit.cellID(), (long long)ahit.cellID(),
                        (long long)ahit.energyDeposit() * 100, 0);
        rawhits->push_back(rawhit);
      }
      return StatusCode::SUCCESS;
    }

    DataHandle<dd4pod::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::RawCalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
  };
  DECLARE_COMPONENT(ExampleCaloDigi)

  // class DataProducerProp : public GaudiAlgorithm {
  // public:
  //  //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
  //  DataProducerProp(const std::string& name, ISvcLocator* svcLoc)
  //      : GaudiAlgorithm(name, svcLoc) {}
  //  StatusCode initialize() override {
  //    f_counter = m_starting_value.value();
  //    return StatusCode::SUCCESS;
  //  }
  //  StatusCode execute() override {
  //    m_vec.put(ThreeVecEx{f_counter, f_counter + 1, f_counter + 2});

  //    info() << "executing DataProducer (prop): " << f_counter << " ..."
  //           << endmsg;
  //    f_counter++;
  //    return StatusCode::SUCCESS;
  //  }

  // private:
  //  Gaudi::Property<int>      m_starting_value{this, "StartingValue", 0};
  //  AnyDataHandle<ThreeVecEx> m_vec{"/Event/UnknownVec2",
  //                                  Gaudi::DataHandle::Writer, this};
  //  int                       f_counter = 0;
  //};

  // DECLARE_COMPONENT(DataProducerProp)

  // class DataConsumerProp : public GaudiAlgorithm {
  // public:
  //  // using GaudiAlgorithm::GaudiAlgorithm;
  //  DataConsumerProp(const std::string& name, ISvcLocator* svcLoc)
  //      : GaudiAlgorithm(name, svcLoc) {}
  //  StatusCode execute() override {
  //    info() << "executing DataConsumer (prop): {" << *m_vec.get() << "}"
  //           << endmsg;
  //    return StatusCode::SUCCESS;
  //  }

  // private:
  //  AnyDataHandle<ThreeVecEx> m_vec{"/Event/UnknownVec2",
  //                                  Gaudi::DataHandle::Reader, this};
  //};

  // DECLARE_COMPONENT(DataConsumerProp)

  } // namespace Examples
} // namespace Gaudi

