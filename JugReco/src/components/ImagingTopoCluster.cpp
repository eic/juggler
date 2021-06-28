/*
 *  Topological Cell Clustering Algorithm for Imaging Calorimetry
 *  1. group all the adjacent pixels
 *
 *  Author: Chao Peng (ANL), 06/02/2021
 *  References: https://arxiv.org/pdf/1603.02934.pdf
 *
 */
#include <algorithm>
#include "fmt/format.h"

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
#include "DD4hep/BitFieldCoder.h"

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/ImagingPixelCollection.h"
#include "eicd/ImagingLayerCollection.h"
#include "eicd/ImagingClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

class ImagingTopoCluster : public GaudiAlgorithm {
public:
    // maximum difference in layer numbers that can be considered as neighbours
    Gaudi::Property<int> m_adjLayerDiff{this, "adjLayerDiff", 1};
    // maximum distance of local (x, y) to be considered as neighbors at the same layer
    Gaudi::Property<std::vector<double>> u_localRanges{this, "localRanges", {1.0*mm, 1.0*mm}};
    // maximum distance of global (eta, phi) to be considered as neighbors at different layers
    Gaudi::Property<std::vector<double>> u_adjLayerRanges{this, "adjLayerRanges", {0.01*M_PI, 0.01*M_PI}};
    // maximum global distance to be considered as neighbors in different sectors
    Gaudi::Property<double> m_adjSectorDist{this, "adjSectorDist", 1.0*cm};
    // minimum cluster center energy (to be considered as a seed for cluster)
    Gaudi::Property<double> m_minClusterCenterEdep{this, "minClusterCenterEdep", 0.};
    // minimum cluster energy (to save this cluster)
    Gaudi::Property<double> m_minEdep{this, "minClusterEdep", 0.5*MeV};
    // minimum number of hits (to save this cluster)
    Gaudi::Property<int> m_minNhits{this, "minClusterNhits", 10};
    // input hits collection
    DataHandle<eic::ImagingPixelCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    // output layer collection
    DataHandle<eic::ImagingLayerCollection>
        m_outputLayerCollection{"outputLayerCollection", Gaudi::DataHandle::Writer, this};
    // output cluster collection
    DataHandle<eic::ImagingClusterCollection>
        m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer, this};

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ImagingTopoCluster(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputHitCollection",       m_inputHitCollection,       "");
        declareProperty("outputClusterCollection",  m_outputClusterCollection,  "");
        declareProperty("outputLayerCollection",  m_outputClusterCollection,  "");
    }

    StatusCode initialize() override
    {
        if (GaudiAlgorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }

        // check local clustering range
        if (u_localRanges.size() < 2) {
            error() << "Need 2-dimensional ranges for same-layer clustering, only have "
                    << u_localRanges.size() << endmsg;
            return StatusCode::FAILURE;
        }
        info() << "Same layer hits group ranges"
               << " (" << u_localRanges[0]/mm << " mm, " << u_localRanges[1]/mm << " mm)"
               << endmsg;

        // check adjacent layer clustering range
        if (u_adjLayerRanges.size() < 2) {
            error() << "Need 2-dimensional ranges for adjacent-layer clustering, only have "
                    << u_adjLayerRanges.size() << endmsg;
            return StatusCode::FAILURE;
        }
        info() << "Same layer hits group ranges"
               << " (" << u_adjLayerRanges[0]/M_PI*1000. << " mrad, "
               << u_adjLayerRanges[1]/M_PI*1000. << " mrad)"
               << endmsg;

        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // input collections
	    const auto &hits = *m_inputHitCollection.get();
        // Create output collections
        auto &clusters = *m_outputClusterCollection.createAndPut();
        auto &layers = *m_outputLayerCollection.createAndPut();

        // group neighboring hits
        std::vector<bool> visits(hits.size(), false);
        std::vector<std::vector<eic::ImagingPixel>> groups;
        for (size_t i = 0; i < hits.size(); ++i)
        {
            // already in a group, or not energetic enough to form a cluster
            if (visits[i] || hits[i].edep() < m_minClusterCenterEdep) {
                continue;
            }
            groups.emplace_back();
            // create a new group, and group all the neighboring hits
            dfs_group(groups.back(), i, hits, visits);
        }
        debug() << "we have " << groups.size() << " groups of hits" << endmsg;

        for (auto &group : groups) {
            add_cluster_layers(group, clusters, layers);
        }

        return StatusCode::SUCCESS;
    }

private:
    template<typename T> static inline T pow2(const T &x) { return x*x; }

    // helper function to group hits
    bool is_neighbor(const eic::ConstImagingPixel &h1, const eic::ConstImagingPixel &h2) const
    {
        // different sectors, simple distance check
        if (h1.sectorID() != h2.sectorID()) {
            return std::sqrt(pow2(h1.x() - h2.x()) + pow2(h1.y() - h2.y()) + pow2(h1.z() - h2.z())) <= m_adjSectorDist;
        }

        // layer check
        int ldiff = std::abs(h1.layerID() - h2.layerID());
        // same layer, check local positions
        if (!ldiff) {
            return (std::abs(h1.local_x() - h2.local_x()) <= u_localRanges[0]) &&
                   (std::abs(h1.local_y() - h2.local_y()) <= u_localRanges[1]);
        } else if (ldiff <= m_adjLayerDiff) {
            return (std::abs(h1.eta() - h2.eta()) <= u_adjLayerRanges[0]) &&
                   (std::abs(h1.phi() - h2.phi()) <= u_adjLayerRanges[1]);
        }

        // not in adjacent layers
        return false;
    }

    // grouping function with Depth-First Search
    void dfs_group(std::vector<eic::ImagingPixel> &group, int idx,
                   const eic::ImagingPixelCollection &hits, std::vector<bool> &visits) const
    {
        auto hit = hits[idx];
        group.push_back(hit);
        visits[idx] = true;
        for(size_t i = 0; i < hits.size(); ++i)
        {
            // visited, or not a neighbor
            if(visits[i] || !is_neighbor(hit, hits[i])) {
                continue;
            }
            dfs_group(group, i, hits, visits);
        }
    }

    void add_cluster_layers(std::vector<eic::ImagingPixel> &group, eic::ImagingClusterCollection &clusters,
                            eic::ImagingLayerCollection &layers)
    {
        // criteria
        if ((int) group.size() < m_minNhits) {
            return;
        }
        double edep = 0.;
        for (auto &hit : group) {
            edep += hit.edep();
        }
        if (edep < m_minEdep) {
            return;
        }

        // build cluster
        auto cl = clusters.create();
        int cid = clusters.size();
        cl.edep(edep);
        cl.nhits(group.size());
        // group hits by layers
        std::map<int, std::vector<size_t>> layer_group;
        for (size_t i = 0; i < group.size(); ++i) {
            auto &hit = group[i];
            int hid = cl.hits_size();
            cl.addhits(hit);
            hit.clusterID(cid);
            hit.hitID(hid);
            auto it = layer_group.find(hit.layerID());
            if (it != layer_group.end()) {
                it->second.push_back(i);
            } else {
                layer_group[hit.layerID()] = {i};
            }
        }

        for (auto &it : layer_group) {
            auto ly = layers.create();
            ly.clusterID(cid);
            ly.layerID(it.first);
            ly.nhits(it.second.size());
            double ledep = 0.;
            for (auto &i : it.second) {
                ly.addhits(group[i]);
                ledep += group[i].edep();
            }
            ly.edep(ledep);
        }
    }

};

DECLARE_COMPONENT(ImagingTopoCluster)

} // namespace Jug::Reco

