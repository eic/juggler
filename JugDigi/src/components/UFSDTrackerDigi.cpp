#include <algorithm>
#include <cmath>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"

// Event Model related classes
// edm4hep's tracker hit is the input collectiopn
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimTrackerHitCollection.h"
// eicd's RawTrackerHit is the output
#include "eicd/RawTrackerHitCollection.h"

namespace Jug::Digi {

  /** UFSD detector digitization.
   *
   * \ingroup digi
   */
  class UFSDTrackerDigi : public GaudiAlgorithm {
  public:
    Gaudi::Property<double>                  m_timeResolution{this, "timeResolution", 10.};
    Gaudi::Property<double>                  m_threshold{this, "threshold", 0. * Gaudi::Units::keV};
    Rndm::Numbers                            m_gaussDist;
    DataHandle<edm4hep::SimTrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                  this};
    DataHandle<eicd::RawTrackerHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                   this};

  public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    UFSDTrackerDigi(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
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
      auto simhits = m_inputHitCollection.get();
      // Create output collections
      auto rawhits = m_outputHitCollection.createAndPut();
      // eicd::RawTrackerHitCollection* rawHitCollection = new eicd::RawTrackerHitCollection();
      std::map<long long, int> cell_hit_map;
      for (const auto& ahit : *simhits) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "--------------------" << ahit.getCellID() << endmsg;
          debug() << "Hit in cellID = " << ahit.getCellID() << endmsg;
          debug() << "     position = (" << ahit.getPosition().x  << "," << ahit.getPosition().y <<","<< ahit.getPosition().z << ")" << endmsg;
          debug() << "    xy_radius = " << std::hypot(ahit.getPosition().x  , ahit.getPosition().y ) << endmsg;
          debug() << "     momentum = (" << ahit.getMomentum().x  << "," << ahit.getMomentum().y <<","<< ahit.getMomentum().z << ")" << endmsg;
        }
        if (ahit.getEDep() * Gaudi::Units::keV < m_threshold) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "         edep = " << ahit.getEDep() << " (below threshold of " << m_threshold / Gaudi::Units::keV << " keV)" << endmsg;
          }          
          continue;
        } else {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "         edep = " << ahit.getEDep() << endmsg;
          }
        }
        // std::cout << ahit << "\n";
        if (cell_hit_map.count(ahit.getCellID()) == 0) {
          cell_hit_map[ahit.getCellID()] = rawhits->size();
          eicd::RawTrackerHit rawhit(0 /* TBDeleted */,
                                    ahit.getCellID(),
                                    ahit.getMCParticle().getTime() * 1e6 + m_gaussDist() * 1e3, // ns->fs
                                    std::llround(ahit.getEDep() * 1e6));
          rawhits->push_back(rawhit);
        } else {
          auto hit = (*rawhits)[cell_hit_map[ahit.getCellID()]];
          hit.time(ahit.getMCParticle().getTime() * 1e6 + m_gaussDist() * 1e3);
          auto ch = hit.charge();
          hit.charge(ch + std::llround(ahit.getEDep() * 1e6));
        }
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(UFSDTrackerDigi)

} // namespace Jug::Digi
