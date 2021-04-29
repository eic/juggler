/*
 *  Topological Cell Clustering Algorithm for Sampling Calorimetry
 *  1. group all the adjacent modules
 *
 *  Author: Chao Peng (ANL), 04/06/2021
 *  References: https://arxiv.org/pdf/1603.02934.pdf
 *
 */
#include <algorithm>

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
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  class TopologicalCellCluster : public GaudiAlgorithm {
  public:
    Gaudi::Property<std::vector<double>> u_groupRanges{this, "groupRanges", {1.0, 1.0, 1.0}};
    Gaudi::Property<double> m_minClusterCenterEdep{this, "minClusterCenterEdep", 0.5*MeV};
    DataHandle<eic::CalorimeterHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ClusterCollection>
        m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer, this};
    DataHandle<eic::CalorimeterHitCollection>
        m_splitHitCollection{"splitHitCollection", Gaudi::DataHandle::Writer, this};

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    TopologicalCellCluster(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputHitCollection",       m_inputHitCollection,       "");
        declareProperty("outputClusterCollection",  m_outputClusterCollection,  "");
        declareProperty("splitHitCollection",       m_splitHitCollection,       "");
    }

    StatusCode initialize() override
    {
        if (GaudiAlgorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }

        if (u_groupRanges.size() < 3) {
            error() << "Need 3-dimensional ranges for clustering, only have " << u_groupRanges.size() << endmsg;
            return StatusCode::FAILURE;
        }
        info() << "Cluster group ranges are"
               << " (" << u_groupRanges[0] << ", " << u_groupRanges[1] << ", " << u_groupRanges[2] << ")"
               << endmsg;
        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // input collections
	    const auto &hits = *m_inputHitCollection.get();
        // Create output collections
        auto &clusters = *m_outputClusterCollection.createAndPut();
        auto &split_hits = *m_splitHitCollection.createAndPut();

        // group neighboring hits
        std::vector<bool> visits(hits.size(), false);
        std::vector<std::vector<eic::CalorimeterHit>> groups;
        for (size_t i = 0; i < hits.size(); ++i)
        {
            // already in a group
            if (visits[i]) {
                continue;
            }
            groups.emplace_back();
            // create a new group, and group all the neighboring hits
            dfs_group(groups.back(), i, hits, visits);
        }
        debug() << "we have " << groups.size() << " groups of hits" << endmsg;

        // TODO: add splitting
        for (auto &group : groups) {
            auto cl = clusters.create();
            for (auto &hit : group) {
                cl.addhits(hit);
            }
        }

        return StatusCode::SUCCESS;
    }

private:
    // helper function to group hits
    inline bool is_neighbor(const eic::ConstCalorimeterHit &h1, const eic::ConstCalorimeterHit &h2) const
    {
        return (std::abs(h1.x() - h2.x()) <= u_groupRanges[0]) &&
               (std::abs(h1.y() - h2.y()) <= u_groupRanges[1]) &&
               (std::abs(h1.z() - h2.z()) <= u_groupRanges[2]);
    }

    // grouping function with Depth-First Search
    void dfs_group(std::vector<eic::CalorimeterHit> &group, int idx,
                   const eic::CalorimeterHitCollection &hits, std::vector<bool> &visits) const
    {
        auto hit = hits[idx];
        group.push_back(hit);
        visits[idx] = true;
        for(size_t i = 0; i < hits.size(); ++i)
        {
            if(visits[i] || !is_neighbor(hit, hits[i])) {
                continue;
            }
            dfs_group(group, i, hits, visits);
        }
    }

  };

DECLARE_COMPONENT(TopologicalCellCluster)

} // namespace Jug::Reco

