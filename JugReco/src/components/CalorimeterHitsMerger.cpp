// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Jihee Kim, Sylvester Joosten, Whitney Armstrong, Wouter Deconinck

/*
 *  An algorithm to group readout hits from a calorimeter
 *  Energy is summed
 *
 *  Author: Chao Peng (ANL), 03/31/2021
 */
#include <algorithm>
#include <bitset>
#include <tuple>
#include <unordered_map>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "fmt/format.h"
#include "fmt/ranges.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

/** Calorimeter hits merging algorithm.
 *
 *  An algorithm to group readout hits from a calorimeter
 *  Energy is summed
 *
 *  \ingroup reco
 */
class CalorimeterHitsMerger : public GaudiAlgorithm {
private:
  Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
  Gaudi::Property<std::string> m_readout{this, "readoutClass", ""};
  // field names to generate id mask, the hits will be grouped by masking the field
  Gaudi::Property<std::vector<std::string>> u_fields{this, "fields", {"layer"}};
  // reference field numbers to locate position for each merged hits group
  Gaudi::Property<std::vector<int>> u_refs{this, "fieldRefNumbers", {}};
  DataHandle<eicd::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::CalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                  this};

  SmartIF<IGeoSvc> m_geoSvc;
  uint64_t id_mask{0}, ref_mask{0};

public:
  CalorimeterHitsMerger(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputHitCollection", m_outputHitCollection, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    m_geoSvc = service(m_geoSvcName);
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }

    if (m_readout.value().empty()) {
      error() << "readoutClass is not provided, it is needed to know the fields in readout ids" << endmsg;
      return StatusCode::FAILURE;
    }

    try {
      auto id_desc = m_geoSvc->detector()->readout(m_readout).idSpec();
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
      return StatusCode::FAILURE;
    }
    id_mask = ~id_mask;
    info() << fmt::format("ID mask in {:s}: {:#064b}", m_readout, id_mask) << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& inputs = *m_inputHitCollection.get();
    // Create output collections
    auto& outputs = *m_outputHitCollection.createAndPut();

    // find the hits that belong to the same group (for merging)
    std::unordered_map<long long, std::vector<eicd::CalorimeterHit>> merge_map;
    for (const auto& h : inputs) {
      int64_t id = h.getCellID() & id_mask;
      // use the reference field position
      auto it = merge_map.find(id);
      if (it == merge_map.end()) {
        merge_map[id] = {h};
      } else {
        it->second.push_back(h);
      }
    }

    // sort hits by energy from large to small
    std::for_each(merge_map.begin(), merge_map.end(), [](auto& it) {
      std::sort(it.second.begin(), it.second.end(), [](const auto& h1, const auto& h2) {
        return h1.getEnergy() > h2.getEnergy();
      });
    });

    // reconstruct info for merged hits
    // dd4hep decoders
    auto poscon = m_geoSvc->cellIDPositionConverter();
    auto volman = m_geoSvc->detector()->volumeManager();

    for (auto& [id, hits] : merge_map) {
      // reference fields id
      const uint64_t ref_id = id | ref_mask;
      // global positions
      const auto gpos = poscon->position(ref_id);
      // local positions
      auto alignment = volman.lookupDetElement(ref_id).nominal();
      const auto pos = alignment.worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
      debug() << volman.lookupDetElement(ref_id).path() << ", " << volman.lookupDetector(ref_id).path() << endmsg;
      // sum energy
      float energy      = 0.;
      float energyError = 0.;
      float time        = 0;
      float timeError   = 0;
      for (auto& hit : hits) {
        energy += hit.getEnergy();
        energyError += hit.getEnergyError() * hit.getEnergyError();
        time += hit.getTime();
        timeError += hit.getTimeError() * hit.getTimeError();
      }
      energyError = sqrt(energyError);
      time /= hits.size();
      timeError = sqrt(timeError) / hits.size();

      const auto& href = hits.front();

      // create const vectors for passing to hit initializer list
      const decltype(eicd::CalorimeterHitData::position) position(
        gpos.x() / dd4hep::mm, gpos.y() / dd4hep::mm, gpos.z() / dd4hep::mm
      );
      const decltype(eicd::CalorimeterHitData::local) local(
        pos.x(), pos.y(), pos.z()
      );

      outputs.push_back(
          eicd::CalorimeterHit{href.getCellID(),
                              energy,
                              energyError,
                              time,
                              timeError,
                              position,
                              href.getDimension(),
                              href.getSector(),
                              href.getLayer(),
                              local}); // Can do better here? Right now position is mapped on the central hit
    }

    if (msgLevel(MSG::DEBUG)) {
      debug() << "Size before = " << inputs.size() << ", after = " << outputs.size() << endmsg;
    }

    return StatusCode::SUCCESS;
  }

}; // class CalorimeterHitsMerger

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(CalorimeterHitsMerger)

} // namespace Jug::Reco
