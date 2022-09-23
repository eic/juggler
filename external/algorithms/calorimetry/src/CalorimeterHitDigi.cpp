// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Wouter Deconinck, Sylvester Joosten

// A general digitization for CalorimeterHit from simulation
// 1. Smear energy deposit with a/sqrt(E/GeV) + b + c/E or a/sqrt(E/GeV) (relative value)
// 2. Digitize the energy with dynamic ADC range and add pedestal (mean +- sigma)
// 3. Time conversion with smearing resolution (absolute value)
// 4. Signal is summed if the SumFields are provided
//
// Author: Chao Peng
// Date: 06/02/2021

#include <algorithms/calorimetry/CalorimeterHitDigi.h>

#include <algorithm>
#include <cmath>
#include <unordered_map>

#include "fmt/format.h"
#include "fmt/ranges.h"

// Event Model related classes
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4eic/RawCalorimeterHitCollection.h"
#include "edm4eic/RawCalorimeterHitData.h"

namespace algorithms::calorimetry {

void CalorimeterHitDigi::init() {
  // set energy resolution numbers
  for (size_t i = 0; i < u_eRes.size() && i < 3; ++i) {
    eRes[i] = u_eRes[i];
  }

  // need signal sum
  if (!u_fields.value().empty()) {
    // sanity checks
    if (m_readout.value().empty()) {
      error() << "readoutClass is not provided, it is needed to know the fields in readout ids"
              << endmsg;
      return;
    }

    // get decoders
    try {
      auto id_desc = m_geoSvc.detector()->readout(m_readout).idSpec();
      id_mask      = 0;
      std::vector<std::pair<std::string, int>> ref_fields;
      for (size_t i = 0; i < u_fields.size(); ++i) {
        id_mask |= id_desc.field(u_fields[i])->mask();
        // use the provided id number to find ref cell, or use 0
        int ref = i < u_refs.size() ? u_refs[i] : 0;
        ref_fields.emplace_back(u_fields[i], ref);
      }
      ref_mask = id_desc.encode(ref_fields);
      // debug() << fmt::format("Referece id mask for the fields {:#064b}", ref_mask) << endmsg;
    } catch (...) {
      error() << "Failed to load ID decoder for " << m_readout << endmsg;
      return;
    }
    id_mask = ~id_mask;
    info() << fmt::format("ID mask in {:s}: {:#064b}", m_readout, id_mask) << endmsg;
  }
}

void CalorimeterHitDigi::process(
    const CalorimeterHitDigi::Input& input,
    const CalorimeterHitDigi::Output& output
) {
  if (!u_fields.value().empty()) {
    signal_sum_digi(input, output);
  } else {
    single_hits_digi(input, output);
  }
}

void CalorimeterHitDigi::single_hits_digi(
    const CalorimeterHitDigi::Input& input,
    const CalorimeterHitDigi::Output& output
) {
  const auto [simhits] = input;
  auto [rawhits]       = output;

  for (const auto& ahit : *simhits) {
    // Note: juggler internal unit of energy is GeV
    const double eDep    = ahit.getEnergy();

    // apply additional calorimeter noise to corrected energy deposit
    const double eResRel = (eDep > m_threshold)
        ? m_randomSvc.normal() * std::sqrt(
              std::pow(eRes[0] / std::sqrt(eDep), 2) +
              std::pow(eRes[1], 2) +
              std::pow(eRes[2] / (eDep), 2)
          )
        : 0;

    const double ped    = m_pedMeanADC + m_randomSvc.normal() * m_pedSigmaADC;
    const long long adc = std::llround(ped +  eDep * (m_corrMeanScale + eResRel) / m_dyRangeADC * m_capADC);

    double time = std::numeric_limits<double>::max();
    for (const auto& c : ahit.getContributions()) {
      if (c.getTime() <= time) {
        time = c.getTime();
      }
    }
    const long long tdc = std::llround((time + m_randomSvc.normal() * m_tRes) / m_resolutionTDC);

    edm4eic::RawCalorimeterHit rawhit(
      ahit.getCellID(),
      (adc > m_capADC.value() ? m_capADC.value() : adc),
      tdc
    );
    rawhits->push_back(rawhit);
  }
}

void CalorimeterHitDigi::signal_sum_digi(
    const CalorimeterHitDigi::Input& input,
    const CalorimeterHitDigi::Output& output
) {
  const auto [simhits] = input;
  auto [rawhits]       = output;

  // find the hits that belong to the same group (for merging)
  std::unordered_map<long long, std::vector<edm4hep::SimCalorimeterHit>> merge_map;
  for (const auto &ahit : *simhits) {
    int64_t hid = (ahit.getCellID() & id_mask) | ref_mask;
    auto    it  = merge_map.find(hid);

    if (it == merge_map.end()) {
      merge_map[hid] = {ahit};
    } else {
      it->second.push_back(ahit);
    }
  }

  // signal sum
  for (auto &[id, hits] : merge_map) {
    double edep     = hits[0].getEnergy();
    double time     = hits[0].getContributions(0).getTime();
    double max_edep = hits[0].getEnergy();
    // sum energy, take time from the most energetic hit
    for (size_t i = 1; i < hits.size(); ++i) {
      edep += hits[i].getEnergy();
      if (hits[i].getEnergy() > max_edep) {
        max_edep = hits[i].getEnergy();
        for (const auto& c : hits[i].getContributions()) {
          if (c.getTime() <= time) {
            time = c.getTime();
          }
        }
      }
    }

    // safety check
    const double eResRel = (edep > m_threshold)
        ? m_randomSvc.normal() * eRes[0] / std::sqrt(edep) +
          m_randomSvc.normal() * eRes[1] +
          m_randomSvc.normal() * eRes[2] / edep
        : 0;

    double    ped     = m_pedMeanADC + m_randomSvc.normal() * m_pedSigmaADC;
    unsigned long long adc     = std::llround(ped + edep * (1. + eResRel) / m_dyRangeADC * m_capADC);
    unsigned long long tdc     = std::llround((time + m_randomSvc.normal() * m_tRes) / m_resolutionTDC);

    edm4eic::RawCalorimeterHit rawhit(
      id,
      (adc > m_capADC.value() ? m_capADC.value() : adc),
      tdc
    );
    rawhits->push_back(rawhit);
  }
}

} // namespace algorithms::calorimetry
