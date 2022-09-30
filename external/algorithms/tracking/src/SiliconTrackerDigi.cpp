// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten

#include <algorithms/tracking/SiliconTrackerDigi.h>

#include <algorithm>
#include <cmath>

namespace algorithms::tracking {

void SiliconTrackerDigi::init() {}

void SiliconTrackerDigi::process(const SiliconTrackerDigi::Input& input,
                                 const SiliconTrackerDigi::Output& output) {
  const auto [simhits] = input;
  auto [rawhits]       = output;

  std::map<long long, int> cell_hit_map;
  for (const auto& ahit : *simhits) {
    debug() << "--------------------" << ahit.getCellID() << endmsg;
    debug() << "Hit in cellID = " << ahit.getCellID() << endmsg;
    debug() << "     position = (" << ahit.getPosition().x << "," << ahit.getPosition().y << ","
            << ahit.getPosition().z << ")" << endmsg;
    debug() << "    xy_radius = " << std::hypot(ahit.getPosition().x, ahit.getPosition().y)
            << endmsg;
    debug() << "     momentum = (" << ahit.getMomentum().x << "," << ahit.getMomentum().y << ","
            << ahit.getMomentum().z << ")" << endmsg;
    if (ahit.getEDep() * dd4hep::keV < m_threshold) {
      debug() << "         edep = " << ahit.getEDep() << " (below threshold of "
              << m_threshold / dd4hep::keV << " keV)" << endmsg;
      continue;
    } else {
      debug() << "         edep = " << ahit.getEDep() << endmsg;
    }
    if (cell_hit_map.count(ahit.getCellID()) == 0) {
      cell_hit_map[ahit.getCellID()] = rawhits->size();
      edm4eic::RawTrackerHit rawhit(ahit.getCellID(),
                                    ahit.getMCParticle().getTime() * 1e6 +
                                        m_rng.gaussian() * 1e3, // ns->fs
                                    std::llround(ahit.getEDep() * 1e6));
      rawhits->push_back(rawhit);
    } else {
      auto hit = (*rawhits)[cell_hit_map[ahit.getCellID()]];
      hit.setTimeStamp(ahit.getMCParticle().getTime() * 1e6 + m_rng.gaussian() * 1e3);
      auto ch = hit.getCharge();
      hit.setCharge(ch + std::llround(ahit.getEDep() * 1e6));
    }
  }
}

} // namespace algorithms::tracking
