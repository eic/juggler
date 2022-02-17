#include <algorithm>
#include <cmath>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
// edm4hep's tracker hit is the input collectiopn
#include "edm4hep/TrackerHitCollection.h"
// eicd's RawTrackerHit is the output
#include "eicd/RawTrackerHitCollection.h"

namespace Jug::Digi {

  /** Silicon detector digitization.
   *
   * \ingroup digi
   */
  class SiliconTrackerDigi : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:
    Gaudi::Property<double>                  m_timeResolution{this, "timeResolution", 10}; // todo : add units
    Gaudi::Property<double>                  m_threshold{this, "threshold", 0. * Gaudi::Units::keV};
    Rndm::Numbers                            m_gaussDist;
    DataHandle<edm4hep::TrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                  this};
    DataHandle<eic::RawTrackerHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                   this};

  public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    SiliconTrackerDigi(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin(name, info())
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
      const edm4hep::TrackerHitCollection* simhits = m_inputHitCollection.get();
      // Create output collections
      auto rawhits = m_outputHitCollection.createAndPut();
      // eic::RawTrackerHitCollection* rawHitCollection = new eic::RawTrackerHitCollection();
      std::map<long long, int> cell_hit_map;
      int ID = 0;
      for (const auto& ahit : *simhits) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "--------------------" << ahit.cellID() << endmsg;
          debug() << "Hit in cellID = " << ahit.cellID() << endmsg;
          debug() << "     position = (" << ahit.position().x  << "," << ahit.position().y <<","<< ahit.position().z << ")" << endmsg;
          debug() << "    xy_radius = " << std::hypot(ahit.position().x  , ahit.position().y ) << endmsg;
          debug() << "     momentum = (" << ahit.momentum().x  << "," << ahit.momentum().y <<","<< ahit.momentum().z << ")" << endmsg;
        }
        if (ahit.energyDeposit() * Gaudi::Units::keV < m_threshold) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "         edep = " << ahit.energyDeposit() << " (below threshold of " << m_threshold / Gaudi::Units::keV << " keV)" << endmsg;
          }          
          continue;
        } else {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "         edep = " << ahit.energyDeposit() << endmsg;
          }
        }
        // std::cout << ahit << "\n";
        if (cell_hit_map.count(ahit.cellID()) == 0) {
          cell_hit_map[ahit.cellID()] = rawhits->size();
          eic::RawTrackerHit rawhit({ID++, algorithmID()},
                                    (int64_t)ahit.cellID(),
                                    ahit.truth().time * 1e6 + m_gaussDist() * 1e3, // ns->fs
                                    std::llround(ahit.energyDeposit() * 1e6));
          rawhits->push_back(rawhit);
        } else {
          auto hit = (*rawhits)[cell_hit_map[ahit.cellID()]];
          hit.time(ahit.truth().time * 1e6 + m_gaussDist() * 1e3);
          auto ch = hit.charge();
          hit.charge(ch + std::llround(ahit.energyDeposit() * 1e6));
        }
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(SiliconTrackerDigi)

} // namespace Jug::Digi
