// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Sylvester Joosten, Wouter Deconinck, Chao, Whitney Armstrong

// Reconstruct digitized outputs, paired with Jug::Digi::CalorimeterHitDigi
// Author: Chao Peng
// Date: 06/14/2021

#include <algorithms/calorimetry/CalorimeterHitReco.h>

#include <algorithm>
#include <functional>
#include <map>

#include <fmt/format.h>
#include <fmt/ranges.h>

namespace algorithms::calorimetry {

void CalorimeterHitReco::init() {

  // unitless conversion
  dyRangeADC = m_dyRangeADC.value() / dd4hep::GeV;

  // threshold for firing
  thresholdADC = m_thresholdFactor.value() * m_pedSigmaADC.value() + m_thresholdValue.value();

  // TDC channels to timing conversion
  stepTDC = dd4hep::ns / m_resolutionTDC.value();

  // do not get the layer/sector ID if no readout class provided
  if (m_readout.value().empty()) {
    return;
  }

  auto id_spec = m_geo->detector()->readout(m_readout).idSpec();
  try {
    id_dec = id_spec.decoder();
    if (!m_sectorField.value().empty()) {
      sector_idx = id_dec->index(m_sectorField);
      info() << "Find sector field " << m_sectorField.value() << ", index = " << sector_idx << endmsg;
    }
    if (!m_layerField.value().empty()) {
      layer_idx = id_dec->index(m_layerField);
      info() << "Find layer field " << m_layerField.value() << ", index = " << sector_idx << endmsg;
    }
  } catch (...) {
    error() << "Failed to load ID decoder for " << m_readout << endmsg;
    return;
  }

  // local detector name has higher priority
  if (!m_localDetElement.value().empty()) {
    try {
      local = m_geo->detector()->detector(m_localDetElement.value());
      info() << "Local coordinate system from DetElement " << m_localDetElement.value() << endmsg;
    } catch (...) {
      error() << "Failed to locate local coordinate system from DetElement " << m_localDetElement.value() << endmsg;
      return;
    }
    // or get from fields
  } else {
    std::vector<std::pair<std::string, int>> fields;
    for (auto& f : u_localDetFields.value()) {
      fields.emplace_back(f, 0);
    }
    local_mask = id_spec.get_mask(fields);
    // use all fields if nothing provided
    if (fields.empty()) {
      local_mask = ~0;
    }
    info() << fmt::format("Local DetElement mask {:#064b} from fields [{}]", local_mask, fmt::join(fields, ", "))
            << endmsg;
  }

  return;
}

void CalorimeterHitReco::process(const CalorimeterHitReco::Input& input,
                                 const CalorimeterHitReco::Output& output) const {
  const auto rawhits = input;
  auto hits          = output;

  auto converter = m_geo->cellIDPositionConverter();

  // energy time reconstruction
  for (const auto& rh : rawhits) {

    #pragma GCC diagnostic push
    #pragma GCC diagnostic error "-Wsign-conversion"

    // did not pass the zero-suppression threshold
    if (rh.getAmplitude() < m_pedMeanADC + thresholdADC) {
      continue;
    }

    // convert ADC -> energy
    const float energy =
      (((signed)rh.getAmplitude() - (signed)m_pedMeanADC)) / static_cast<float>(m_capADC.value()) * dyRangeADC / m_sampFrac;
    const float time  = rh.getTimeStamp() / stepTDC;

    #pragma GCC diagnostic pop

    const auto cellID = rh.getCellID();
    const int lid = id_dec != nullptr && !m_layerField.value().empty() ? static_cast<int>(id_dec->get(cellID, layer_idx)) : -1;
    const int sid = id_dec != nullptr && !m_sectorField.value().empty() ? static_cast<int>(id_dec->get(cellID, sector_idx)) : -1;
    // global positions
    const auto gpos = converter->position(cellID);
    // local positions
    if (m_localDetElement.value().empty()) {
      auto volman = m_geo->detector()->volumeManager();
      local       = volman.lookupDetElement(cellID & local_mask);
      }
      const auto pos = local.nominal().worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
      // cell dimension
      std::vector<double> cdim;
      // get segmentation dimensions
      if (converter->findReadout(local).segmentation().type() != "NoSegmentation") {
      cdim = converter->cellDimensions(cellID);
      // get volume dimensions (multiply by two to get fullsize)
      } else {
      // Using bounding box instead of actual solid so the dimensions are always in dim_x, dim_y, dim_z
      cdim = converter->findContext(cellID)->volumePlacement().volume().boundingBox().dimensions();
      std::transform(cdim.begin(), cdim.end(), cdim.begin(),
                      std::bind(std::multiplies<double>(), std::placeholders::_1, 2));
      }

      // create const vectors for passing to hit initializer list
      const decltype(edm4eic::CalorimeterHitData::position) position(
        gpos.x() / m_lUnit, gpos.y() / m_lUnit, gpos.z() / m_lUnit
      );
      const decltype(edm4eic::CalorimeterHitData::dimension) dimension(
        cdim[0] / m_lUnit, cdim[1] / m_lUnit, cdim[2] / m_lUnit
      );
      const decltype(edm4eic::CalorimeterHitData::local) local_position(
        pos.x() / m_lUnit, pos.y() / m_lUnit, pos.z() / m_lUnit
      );

      hits.push_back({
          rh.getCellID(), // cellID
          energy,         // energy
          0,              // @TODO: energy error
          time,           // time
          0,              // time error FIXME should be configurable
          position,       // global pos
          dimension,
          // Local hit info
          sid,
          lid,
          local_position, // local pos
      });
  }
}

} // namespace algorithms::calorimetry
