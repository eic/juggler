/*
 *  Island Clustering Algorithm for Calorimeter Blocks
 *  1. group all the adjacent modules
 *  2. split the groups between their local maxima with the energy deposit above <minClusterCenterEdep>
 *
 *  Author: Chao Peng (ANL), 09/27/2020
 *  References:
 *      https://cds.cern.ch/record/687345/files/note01_034.pdf
 *      https://www.jlab.org/primex/weekly_meetings/primexII/slides_2012_01_20/island_algorithm.pdf
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
class CalorimeterIslandCluster : public GaudiAlgorithm
{
public:
    Gaudi::Property<double> m_groupRange{this, "groupRange", 1.8};
    Gaudi::Property<double> m_minClusterCenterEdep{this, "minClusterCenterEdep", 50.0*MeV};
    DataHandle<eic::CalorimeterHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ClusterCollection>
        m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    CalorimeterIslandCluster(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputHitCollection",      m_inputHitCollection,      "");
        declareProperty("outputClusterCollection", m_outputClusterCollection, "");
    }

    StatusCode initialize() override
    {
        if (GaudiAlgorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }
        m_geoSvc = service("GeoSvc");
        if (!m_geoSvc) {
            error() << "Unable to locate Geometry Service. "
                    << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
            return StatusCode::FAILURE;
        }
        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // input collections
	    const auto &hits = *m_inputHitCollection.get();
        // Create output collections
        auto &clusters = *m_outputClusterCollection.createAndPut();

        // group neighboring hits
        std::vector<bool> visits(hits.size(), false);
        eic::ClusterCollection groups;
        for (size_t i = 0; i < hits.size(); ++i)
        {
            // already in a group
            if (visits[i]) {
                continue;
            }
            // create a new group, and group all the neighboring hits
            dfs_group(groups.create(), i, hits, visits);
        }
        info() << "we have " << groups.size() << " groups of hits" << endmsg;

        for (auto &group : groups) {
            auto maxima = find_local_maxima(group);
            auto split_collection = split_group(group, maxima, clusters);
            info() << "hits in a group: " << group.hits_size() <<  ", "
                   << "local maxima: " << maxima.hits_size() << endmsg;
        }

        return StatusCode::SUCCESS;
    }

private:
    // helper function to group hits
    inline bool is_neighbor(const eic::ConstCalorimeterHit &h1, const eic::ConstCalorimeterHit &h2)
    {
        auto pos1 = h1.position();
        auto pos2 = h2.position();
        auto dim1 = m_geoSvc->cellIDPositionConverter()->cellDimensions(h1.cellID0());
        auto dim2 = m_geoSvc->cellIDPositionConverter()->cellDimensions(h2.cellID0());

        // info() << std::abs(pos1.x - pos2.x) << ", " << (dim1[0] + dim2[0])/2. << ", "
        //        << std::abs(pos1.y - pos2.y) << ", " << (dim1[1] + dim2[1])/2. << endmsg;

        return (std::abs(pos1.x - pos2.x) <= (dim1[0] + dim2[0])/2.*m_groupRange) &&
               (std::abs(pos1.y - pos2.y) <= (dim1[1] + dim2[1])/2.*m_groupRange);
    }

    // grouping function with Depth-First Search
    void dfs_group(eic::Cluster group, int idx, const eic::CalorimeterHitCollection &hits, std::vector<bool> &visits)
    {
        auto hit = hits[idx];
        group.addhits(hit);
        visits[idx] = true;
        for(size_t i = 0; i < hits.size(); ++i)
        {
            if(visits[i] || !is_neighbor(hit, hits[i])) {
                continue;
            }
            dfs_group(group, i, hits, visits);
        }
    }

    // find local maxima that above a certain threshold
    eic::Cluster find_local_maxima(const eic::Cluster &group)
    {
        eic::Cluster maxima;
        for(auto &hit : group.hits())
        {
            // not a qualified center
            if(hit.energy() < m_minClusterCenterEdep) {
                continue;
            }

            bool maximum = true;
            for(auto &hit2 : group.hits())
            {
                if(hit == hit2)
                    continue;

                if(is_neighbor(hit, hit2) && hit2.energy() > hit.energy()) {
                    maximum = false;
                    break;
                }
            }

            if(maximum) {
                maxima.addhits(hit);
            }
        }

        return maxima;
    }

    // helper function
    inline void vec_normalize(std::vector<double> &vals) {
        double total = 0.;
        for (auto &val : vals) { total += val; }
        for (auto &val : vals) { val /= total; }
    }

    // split a group of hits according to the local maxima
    eic::CalorimeterHitCollection split_group(eic::Cluster group, const eic::Cluster &maxima,
                                              eic::ClusterCollection &clusters)
    {
        // to persistify the split hits
        eic::CalorimeterHitCollection scoll;
        // special cases
        if (maxima.hits_size() == 0) {
            return scoll;
        } else if (maxima.hits_size() == 1) {
            clusters.push_back(group.clone());
            return scoll;
        }

        // distance reference
        auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(maxima.hits_begin()->cellID0());
        double dist_ref = dim[0];

        // split between maxima
        std::vector<double> weights(maxima.hits_size());
        std::vector<eic::Cluster> splits(maxima.hits_size());
        size_t i = 0;
        for (auto it = group.hits_begin(); it != group.hits_end(); ++it, ++i) {
            auto hpos = it->position();
            auto hedep = it->energy();
            size_t j = 0;
            // calculate weights for local maxima
            for (auto cit = maxima.hits_begin(); cit != maxima.hits_end(); ++cit, ++j) {
                double energy = cit->energy();
                auto pos = cit->position();
                double dist = std::sqrt(std::pow(pos.x - hpos.x, 2) + std::pow(pos.y - hpos.y, 2));
                weights[j] = std::exp(-dist/dist_ref)*energy;
            }

            // normalize weights
            vec_normalize(weights);

            // ignore small weights
            for (auto &w : weights) {
                if (w < 0.02) { w = 0; }
            }
            vec_normalize(weights);

            // split energy between local maxima
            for (size_t k = 0; k < weights.size(); ++k) {
                double weight = weights[k];
                if (weight <= 1e-6) {
                    continue;
                }

                eic::CalorimeterHit hit(it->cellID0(), it->cellID1(), hedep*weight,
                                        it->time(), it->position(), it->type());
                scoll.push_back(hit);
                splits[k].addhits(hit);
            }
        }

        for (auto &cl : splits) {
            clusters.push_back(cl);
        }
        return scoll;
    }
};

DECLARE_COMPONENT(CalorimeterIslandCluster)

} // namespace Jug::Reco

