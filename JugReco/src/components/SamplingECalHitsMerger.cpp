/*
 *  A clustering algorithm to reduce 3D clustering to 2D (x, y) + 1D (depth)
 *  2D clustering is formed by summing all layers in the same module
 *  1D clustering is formed by summing all modules in the same layer
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

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

class SamplingECalHitsMerger : public GaudiAlgorithm {
public:
    Gaudi::Property<std::vector<std::pair<int, int>>>
        u_cellIDMaskRanges{this, "cellIDMaskRanges", {{0, 31}}};
    DataHandle<eic::CalorimeterHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection>
        m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
    int64_t id_mask;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    SamplingECalHitsMerger(const std::string& name, ISvcLocator* svcLoc)
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

        // build masks from input
        id_mask = 0;
        for (auto &p : u_cellIDMaskRanges) {
            debug() << "masking bit " << p.first << " - " << p.second << endmsg;
            for (int64_t k = p.first; k <= p.second; ++k) {
                id_mask |= (int64_t(1) << k);
            }
        }
        id_mask = ~id_mask;
        debug() << "cellID mask = " << std::bitset<64>(id_mask) << endmsg;

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

        for (auto &h : mhits) {
            debug() << h.cellID() << ": " << h.energy() << endmsg;
        }

        debug() << "Size before = " << hits.size() << ", after = " << mhits.size() << endmsg;

        return StatusCode::SUCCESS;
    }


}; // class SamplingECalHitsMerger

DECLARE_COMPONENT(SamplingECalHitsMerger)

} // namespace Jug::Reco

