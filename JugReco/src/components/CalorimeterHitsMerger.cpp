/*
 *  An algorithm to group readout hits from a calorimeter
 *  Energy is summed
 *
 *  Author: Chao Peng (ANL), 03/31/2021
 */
#include <tuple>
#include <bitset>
#include <algorithm>
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
#include "JugBase/Utilities/UniqueID.hpp"

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
    // Unique identifier for this hit type, based on the algorithm name
    using HitClassificationType = decltype(eic::CalorimeterHitData().type);
    const HitClassificationType m_type;
  public:
    Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
    Gaudi::Property<std::string> m_readout{this, "readoutClass", ""};
    // field names to generate id mask, the hits will be grouped by masking the field
    Gaudi::Property<std::vector<std::string>> u_fields{this, "fields", {"layer"}};
    // reference field numbers to locate position for each merged hits group
    Gaudi::Property<std::vector<int>>         u_refs{this, "fieldRefNumbers", {}};
    DataHandle<eic::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection",
                                                                   Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection> m_outputHitCollection{
        "outputHitCollection", Gaudi::DataHandle::Writer, this};

    SmartIF<IGeoSvc> m_geoSvc;
    uint64_t         id_mask, ref_mask;

    CalorimeterHitsMerger(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc)
      , m_type{uniqueID<HitClassificationType>(name)}
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }

      m_geoSvc = service(m_geoSvcName);
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }

      if (m_readout.value().empty()) {
        error() << "readoutClass is not provided, it is needed to know the fields in readout ids"
                << endmsg;
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
          ref_fields.push_back({u_fields[i], ref});
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

    StatusCode execute() override
    {
      // input collections
      const auto& inputs = *m_inputHitCollection.get();
      // Create output collections
      auto& outputs = *m_outputHitCollection.createAndPut();

      // find the hits that belong to the same group (for merging)
      std::unordered_map<long long, std::vector<eic::ConstCalorimeterHit>> merge_map;
      for (const auto& h : inputs) {
        int64_t id = h.cellID() & id_mask;
        // use the reference field position
        auto it = merge_map.find(id);
        if (it == merge_map.end()) {
          merge_map[id] = {h};
        } else {
          it->second.push_back(h);
        }
      }

      // reconstruct info for merged hits
      // dd4hep decoders
      auto poscon = m_geoSvc->cellIDPositionConverter();
      auto volman = m_geoSvc->detector()->volumeManager();

      int nresults = 0;
      for (auto &[id, hits] : merge_map) {
        // reference fields id
        int64_t ref_id = id | ref_mask;
        // global positions
        auto gpos = poscon->position(ref_id);
        // local positions
        auto alignment = volman.lookupDetElement(ref_id).nominal();
        auto pos = alignment.worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
        // debug() << volman.lookupDetElement(ref_id).path() << ", "
        //         << volman.lookupDetector(ref_id).path() << endmsg;
        // sum energy
        float energy = 0.;
        for (auto &hit : hits) {
          energy += hit.energy();
          // debug() << fmt::format("{:#064b} - {:#064b}, ref: {:#064b}", hit.cellID(), id, ref_id)
          //         << endmsg;
        }
        const auto &href = hits.front();
        outputs.push_back(eic::CalorimeterHit{
                          ref_id,
                          nresults++,
                          href.layer(),
                          href.sector(),
                          m_type,
                          energy,
                          0, //@TODO: energy uncertainty
                          href.time(),
                          {gpos.x() / dd4hep::mm, gpos.y() / dd4hep::mm, gpos.z() / dd4hep::mm},
                          {pos.x() / dd4hep::mm, pos.y() / dd4hep::mm, pos.z() / dd4hep::mm},
                          href.dimension()});
      }

      debug() << "Size before = " << inputs.size() << ", after = " << outputs.size() << endmsg;

      return StatusCode::SUCCESS;
    }

  }; // class CalorimeterHitsMerger

  DECLARE_COMPONENT(CalorimeterHitsMerger)

} // namespace Jug::Reco

