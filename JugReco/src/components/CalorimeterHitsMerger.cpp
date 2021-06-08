/*
 *  An algorithm to group readout hits from a calorimeter
 *  Energy is summed
 *
 *  Author: Chao Peng (ANL), 03/31/2021
 */
#include <bitset>
#include <algorithm>
#include <unordered_map>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/PhysicalConstants.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "fmt/ranges.h"
#include "fmt/format.h"

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

class CalorimeterHitsMerger : public GaudiAlgorithm {
public:
    Gaudi::Property<std::string>                m_geoSvcName{this, "geoServiceName", "GeoSvc"};
    Gaudi::Property<std::string>                m_readout{this, "readoutClass", "EcalBarrelHits"};
    Gaudi::Property<std::vector<std::string>>   m_fields{this, "fields", {"layer"}};
    DataHandle<eic::CalorimeterHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection>
        m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};

    SmartIF<IGeoSvc> m_geoSvc;
    uint64_t id_mask;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    CalorimeterHitsMerger(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputHitCollection",       m_inputHitCollection,       "");
        declareProperty("outputHitCollection",      m_outputHitCollection,    "");
    }

    StatusCode initialize() override
    {
        if (GaudiAlgorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }

        m_geoSvc = service(m_geoSvcName);
        if (!m_geoSvc) {
            error() << "Unable to locate Geometry Service. "
                    << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
            return StatusCode::FAILURE;
        }

        try {
            auto id_desc = m_geoSvc->detector()->readout(m_readout).idSpec();
            id_mask = 0;
            for (auto &f : m_fields) {
                id_mask |= id_desc.field(f)->mask();
            }
        } catch (...) {
            error() << "Failed to load ID decoder for " << m_readout << endmsg;
            return StatusCode::FAILURE;
        }
        id_mask = ~id_mask;
        info() << fmt::format("ID mask for [{:s}] fields in {:s}: {:#064b}",
                              fmt::join(m_fields, ", "), m_readout, id_mask)
               << endmsg;
        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // input collections
	    const auto &hits = *m_inputHitCollection.get();
        // Create output collections
        auto &mhits = *m_outputHitCollection.createAndPut();

        // sum energies that has the same id
        std::unordered_map<long long, size_t> merge_map;
        for (auto &h : hits) {
            auto id = (h.cellID() & id_mask);
            // debug() << h.cellID() << " - " << std::bitset<64>(h.cellID()) << endmsg;
            auto it = merge_map.find(id);
            if (it == merge_map.end()) {
                merge_map[id] = mhits.size();
                mhits.push_back(h.clone());
                debug() << mhits[mhits.size() - 1].cellID() << " - " << std::bitset<64>(id) << endmsg;
            } else {
                mhits[it->second].energy(mhits[it->second].energy() + h.energy());
            }
        }

        /*
        for (auto &h : mhits) {
            debug() << h.cellID() << ": " << h.energy() << endmsg;
        }
        */

        debug() << "Size before = " << hits.size() << ", after = " << mhits.size() << endmsg;

        return StatusCode::SUCCESS;
    }


}; // class CalorimeterHitsMerger

DECLARE_COMPONENT(CalorimeterHitsMerger)

} // namespace Jug::Reco

