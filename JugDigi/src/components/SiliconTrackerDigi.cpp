// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten

#include <algorithm>
#include <cmath>

#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"

#include <k4FWCore/DataHandle.h>

// Event Model related classes
// edm4hep's tracker hit is the input collectiopn
#include "edm4hep/EDM4hepVersion.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
// edm4eic's RawTrackerHit is the output
#include "edm4eic/RawTrackerHitCollection.h"

namespace Jug::Digi {

/** Silicon detector digitization.
 *
 * \ingroup digi
 */
class SiliconTrackerDigi : public Gaudi::Algorithm {
private:
  Gaudi::Property<double> m_timeResolution{this, "timeResolution", 10}; // todo : add units
  Gaudi::Property<double> m_threshold{this, "threshold", 0. * Gaudi::Units::keV};
  Rndm::Numbers m_gaussDist;
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                    this};
  mutable DataHandle<edm4eic::RawTrackerHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                  this};

public:
  //  ill-formed: using Gaudi::Algorithm::GaudiAlgorithm;
  SiliconTrackerDigi(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputHitCollection", m_outputHitCollection, "");
  }
  StatusCode initialize() override {
    if (Gaudi::Algorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    IRndmGenSvc* randSvc = Gaudi::svcLocator()->service<IRndmGenSvc>("RndmGenSvc", true);
    StatusCode sc        = m_gaussDist.initialize(randSvc, Rndm::Gauss(0.0, m_timeResolution.value()));
    if (!sc.isSuccess()) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
  StatusCode execute(const EventContext&) const override {
    // input collection
    const auto* const simhits = m_inputHitCollection.get();
    // Create output collections
    auto* rawhits = m_outputHitCollection.createAndPut();
    // edm4eic::RawTrackerHitCollection* rawHitCollection = new edm4eic::RawTrackerHitCollection();
    std::map<long long, int> cell_hit_map;
    for (const auto& ahit : *simhits) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "--------------------" << ahit.getCellID() << endmsg;
        debug() << "Hit in cellID = " << ahit.getCellID() << endmsg;
        debug() << "     position = (" << ahit.getPosition().x << "," << ahit.getPosition().y << ","
                << ahit.getPosition().z << ")" << endmsg;
        debug() << "    xy_radius = " << std::hypot(ahit.getPosition().x, ahit.getPosition().y) << endmsg;
        debug() << "     momentum = (" << ahit.getMomentum().x << "," << ahit.getMomentum().y << ","
                << ahit.getMomentum().z << ")" << endmsg;
      }
      if (ahit.getEDep() * Gaudi::Units::keV < m_threshold) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "         edep = " << ahit.getEDep() << " (below threshold of " << m_threshold / Gaudi::Units::keV
                  << " keV)" << endmsg;
        }
        continue;
      } else {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "         edep = " << ahit.getEDep() << endmsg;
        }
      }
      if (cell_hit_map.count(ahit.getCellID()) == 0) {
        cell_hit_map[ahit.getCellID()] = rawhits->size();
        rawhits->create(
          ahit.getCellID(),
#if EDM4HEP_BUILD_VERSION >= EDM4HEP_VERSION(0, 99, 0)
          ahit.getParticle().getTime() * 1e6 + m_gaussDist() * 1e3, // ns->fs
#else
          ahit.getMCParticle().getTime() * 1e6 + m_gaussDist() * 1e3, // ns->fs
#endif
          std::llround(ahit.getEDep() * 1e6)
        );
      } else {
        auto hit = (*rawhits)[cell_hit_map[ahit.getCellID()]];
#if EDM4HEP_BUILD_VERSION >= EDM4HEP_VERSION(0, 99, 0)
        hit.setTimeStamp(ahit.getParticle().getTime() * 1e6 + m_gaussDist() * 1e3);
#else
        hit.setTimeStamp(ahit.getMCParticle().getTime() * 1e6 + m_gaussDist() * 1e3);
#endif
        auto ch = hit.getCharge();
        hit.setCharge(ch + std::llround(ahit.getEDep() * 1e6));
      }
    }
    return StatusCode::SUCCESS;
  }
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(SiliconTrackerDigi)

} // namespace Jug::Digi
