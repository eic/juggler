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

  class CalorimeterIslandCluster : public GaudiAlgorithm {
  public:
    Gaudi::Property<bool> m_splitCluster{this, "splitCluster", true};
    Gaudi::Property<double> m_groupRange{this, "groupRange", 1.8};
    Gaudi::Property<double> m_minClusterCenterEdep{this, "minClusterCenterEdep", 50.0*MeV};
    DataHandle<eic::CalorimeterHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ClusterCollection>
        m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer, this};
    DataHandle<eic::CalorimeterHitCollection>
        m_splitHitCollection{"splitHitCollection", Gaudi::DataHandle::Writer, this};

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    CalorimeterIslandCluster(const std::string& name, ISvcLocator* svcLoc)
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

        for (auto &group : groups) {
            auto maxima = find_maxima(group, !m_splitCluster);
            split_group(group, maxima, clusters, split_hits);
            debug() << "hits in a group: " << group.size() <<  ", "
                    << "local maxima: " << maxima.hits_size() << endmsg;
        }

        return StatusCode::SUCCESS;
    }

private:
    // helper function to group hits
    inline bool is_neighbor(const eic::ConstCalorimeterHit &h1, const eic::ConstCalorimeterHit &h2) const
    {
        return (std::abs(h1.local_x() - h2.local_x()) <= (h1.dim_x() + h2.dim_y())/2.*m_groupRange) &&
               (std::abs(h1.local_y() - h2.local_y()) <= (h1.dim_y() + h2.dim_y())/2.*m_groupRange);
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

    // find local maxima that above a certain threshold
    eic::Cluster find_maxima(const std::vector<eic::CalorimeterHit> &group, bool global = false) const
    {
        eic::Cluster maxima;
        if (group.empty()) {
            return maxima;
        }

        if (global) {
            int mpos = 0;
            for (size_t i = 0; i < group.size(); ++i) {
                if (group[mpos].energy() < group[i].energy()) {
                    mpos = i;
                }
            }
            maxima.addhits(group[mpos]);
            return maxima;
        }

        for(auto &hit : group)
        {
            // not a qualified center
            if(hit.energy() < m_minClusterCenterEdep/GeV) {
                continue;
            }

            bool maximum = true;
            for(auto &hit2 : group)
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
    inline void vec_normalize(std::vector<double> &vals) const {
        double total = 0.;
        for (auto &val : vals) { total += val; }
        for (auto &val : vals) { val /= total; }
    }

    // split a group of hits according to the local maxima
    // split_hits is used to persistify the data
    void split_group(const std::vector<eic::CalorimeterHit> &group, const eic::Cluster &maxima,
                     eic::ClusterCollection &clusters, eic::CalorimeterHitCollection &split_hits) const
    {
        // special cases
        if (maxima.hits_size() == 0) {
            return;
        } else if (maxima.hits_size() == 1) {
            auto cl = clusters.create();
            for (auto &hit : group) {
                cl.addhits(hit);
            }
            return;
        }

        // split between maxima
        std::vector<double> weights(maxima.hits_size());
        std::vector<eic::Cluster> splits(maxima.hits_size());
        size_t i = 0;
        for (auto it = group.begin(); it != group.end(); ++it, ++i) {
            auto hedep = it->energy();
            size_t j = 0;
            // calculate weights for local maxima
            for (auto cit = maxima.hits_begin(); cit != maxima.hits_end(); ++cit, ++j) {
                double dist_ref = cit->dim_x();
                double energy = cit->energy();
                double dist = std::sqrt(std::pow(it->local_x() - cit->local_x(), 2)
                                        + std::pow(it->local_y() - cit->local_y(), 2));
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
                auto hit = it->clone();
                hit.energy(hedep*weight);
                hit.type(1);
                split_hits.push_back(hit);
                splits[k].addhits(hit);
            }
        }

        for (auto &cl : splits) {
            clusters.push_back(cl);
        }
        return;
    }
  };

DECLARE_COMPONENT(CalorimeterIslandCluster)

} // namespace Jug::Reco

