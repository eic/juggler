#include <algorithm>

#include "GaudiAlg/Transformer.h"

// Event Model related classes
//#include "GaudiExamples/MyTrack.h"
#include "eicd/RawCalorimeterHitCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"

namespace Jug {
  namespace Digi {

  struct ExampleCaloDigi final
      : Gaudi::Functional::Transformer<eic::RawCalorimeterHitCollection(
            const dd4pod::CalorimeterHitCollection&)> {

    ExampleCaloDigi(const std::string& name, ISvcLocator* pSvc)
        : Transformer(name, pSvc, {KeyValue("InputData", {"FAEC_ShHits"})},
                      KeyValue("OutputData", {"ForwardPreshowerHits"})) {}

    eic::RawCalorimeterHitCollection
    operator()(const dd4pod::CalorimeterHitCollection& in_hits) const override {
      eic::RawCalorimeterHitCollection out_hits;
      for (auto i = in_hits.begin(), end = in_hits.end(); i != end; ++i) {
        out_hits.create((long long)i->cellID(), (long long)i->cellID(),
                        (long long)i->energyDeposit() * 100, 0);
      }
      return out_hits;
    }
    };

    DECLARE_COMPONENT(ExampleCaloDigi)
  } // namespace Examples
} // namespace Gaudi

