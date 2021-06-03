/*
 *  A hits converter to prepare dataset for machine learning
 *  It converts hits with (x, y, z, E) to (E, eta, phi) layer by layer
 *  With a defined grid size and ranges, it merge the hits within one grid and drop-off hits out-of-range
 *  The capacity of each layer is fixed (padding with zeros), and the hits with least energies that exceed the capacity
 *  will be discarded.
 *
 *  Author: Chao Peng (ANL), 05/04/2021
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
#include "JugBase/Utilities/Utils.hpp"

// Event Model related classes
#include "eicd/ImagingPixel.h"
#include "eicd/ImagingPixelCollection.h"
#include "eicd/CalorimeterHitCollection.h"

using namespace Gaudi::Units;

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

namespace Jug::Reco {

class ImagingPixelMerger : public GaudiAlgorithm {
public:
    Gaudi::Property<int> m_nHits{this, "numberOfHits", 20};
    Gaudi::Property<int> m_nLayers{this, "numberOfLayers", 20};
    Gaudi::Property<double> m_etaSize{this, "etaSize", 0.001};
    Gaudi::Property<double> m_phiSize{this, "phiSize", 0.001};
    Gaudi::Property<std::vector<int>> u_layerIDMaskRange{this, "layerIDMaskRange", {}};
    DataHandle<eic::CalorimeterHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ImagingPixelCollection>
        m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ImagingPixelMerger(const std::string& name, ISvcLocator* svcLoc)
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

        // info() << "z_length " << depth << endmsg;
        auto vals = u_layerIDMaskRange.value();
        if (vals.size() != 2) {
            error() << "Need layerIDMaskRange to proceed." << endmsg;
            return StatusCode::FAILURE;
        }

        // build masks from range
        id_shift = vals[0];
        id_mask = 0;
        debug() << "masking bit " << vals[0] << " - " << vals[1] << endmsg;
        for (int64_t k = 0; k <= vals[1] - vals[0]; ++k) {
            id_mask |= (int64_t(1) << k);
        }
        debug() << "layer mask = " << std::bitset<64>(id_mask) << endmsg;

        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // input collections
	    const auto &hits = *m_inputHitCollection.get();
        // Create output collections
        auto &mhits = *m_outputHitCollection.createAndPut();

        // group the hits by layer
        std::vector<std::unordered_map<std::pair<int, int>, double, pair_hash>> group_hits(m_nLayers);
        for (auto &h : hits) {
            auto k = get_subid(h.cellID(), id_mask, id_shift) - 1;
            if ((int) k >= m_nLayers) {
                continue;
            }
            double r = std::sqrt(h.x()*h.x() + h.y()*h.y() + h.z()*h.z());
            double th = std::acos(h.z()/r);
            double eta = -std::log(std::tan(th/2.));
            double phi = std::atan2(h.y(), h.x());
            // debug() << th << ", " << eta << ", " << phi << endmsg;

            auto &layer = group_hits[k];
            auto g = std::pair<int, int>{(eta + 4.)/m_etaSize, (phi + M_PI)/m_phiSize};
            auto it = layer.find(g);
            // merge energy
            if (it != layer.end()) {
                it->second += h.energy();
            } else {
                layer[g] = h.energy();
            }
        }

        // convert to data
        struct GridData { double energy, eta, phi; };
        for (auto [i, layer] : Jug::Utils::Enumerate(group_hits)) {
            std::vector<GridData> grids;
            for (auto &it : layer) {
                double ie = it.first.first;
                double ip = it.first.second;
                grids.push_back(GridData{it.second, ie*m_etaSize - 4., ip*m_phiSize - M_PI});
                // debug() << ie << ", " << ip << ", " << ie*m_etaSize - 4. << ", " << ip*m_phiSize - M_PI << endmsg;
            }
            std::sort(grids.begin(), grids.end(),
                      [] (const GridData &g1, const GridData &g2) { return g1.energy < g2.energy; });

            for (size_t k = 0; k < (size_t) m_nHits; ++k) {
                GridData grid {0., 0., 0.};
                if (k < grids.size()) {
                    grid = grids[k];
                }
                auto h = mhits.create();
                h.edep(grid.energy);
                h.eta(grid.eta);
                h.phi(grid.phi);
                h.layerID(i);
                h.hitID(k);
            }
        }

        return StatusCode::SUCCESS;
    }

private:
    uint64_t id_mask, id_shift;
    // helper function to unfold layer id
    inline uint64_t get_subid(int64_t cid, int64_t mask, int64_t shift) const
    {
        return (cid >> shift) & mask;
    }

}; // class ImagingPixelMerger

DECLARE_COMPONENT(ImagingPixelMerger)

} // namespace Jug::Reco

