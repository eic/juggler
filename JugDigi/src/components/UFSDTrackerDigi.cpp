#include <algorithm>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Property.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
//#include "GaudiExamples/MyTrack.h"
//
// dd4pod's tracker hit is the input collectiopn
#include "dd4pod/TrackerHitCollection.h"
// eicd's RawTrackerHit is the output
#include "eicd/RawTrackerHitCollection.h"

namespace Jug::Digi {

  /** Ultra-fast silicon detector digitization.
   *
   *
   */
  class UFSDTrackerDigi : public GaudiAlgorithm {
  public:
    Gaudi::Property<double>                  m_timeResolution{this, "timeResolution", 10};
    Rndm::Numbers                            m_gaussDist;
    DataHandle<dd4pod::TrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                  this};
    DataHandle<eic::RawTrackerHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                   this};

  public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    UFSDTrackerDigi(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }
    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
      StatusCode   sc      = m_gaussDist.initialize(randSvc, Rndm::Gauss(0.0, m_timeResolution.value()));
      if (!sc.isSuccess()) {
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }
    StatusCode execute() override
    {
      // input collection
      const dd4pod::TrackerHitCollection* simhits = m_inputHitCollection.get();
      // Create output collections
      auto rawhits = m_outputHitCollection.createAndPut();
      // eic::RawTrackerHitCollection* rawHitCollection = new eic::RawTrackerHitCollection();
      std::map<long long, int> cell_hit_map;
      for (const auto& ahit : *simhits) {
        // std::cout << ahit << "\n";
        if (cell_hit_map.count(ahit.cellID()) == 0) {
          cell_hit_map[ahit.cellID()] = rawhits->size();
          eic::RawTrackerHit rawhit((long long)ahit.cellID(),
                                    ahit.truth().time * 1e6 + m_gaussDist() * 1e3, // ns->fs
                                    int(ahit.energyDeposit() * 1e6));
          rawhits->push_back(rawhit);
        } else {
          auto hit = (*rawhits)[cell_hit_map[ahit.cellID()]];
          hit.time(ahit.truth().time * 1e6 + m_gaussDist() * 1e3);
          auto ch = hit.charge();
          hit.charge(ch + int(ahit.energyDeposit() * 1e6));
        }
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(UFSDTrackerDigi)

} // namespace Jug::Digi
