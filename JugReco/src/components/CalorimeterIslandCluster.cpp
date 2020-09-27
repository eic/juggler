/*
 *  Island Clustering Algorithm for Calorimeter Blocks
 *  1. group all the adjacent modules with the energy deposit above <minModuleEdep>
 *  2. split the groups between their local maxima with the energy deposit above <minClusterCenterEdep>
 *  3. reconstruct the clustrers
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
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {
class CalorimeterIslandCluster : public GaudiAlgorithm
{
public:
    Gaudi::Property<double> m_minModuleEdep{this, "minModuleEdep", 0.5*MeV};
    Gaudi::Property<double> m_minClusterCenterEdep{this, "minClusterCenterEdep", 50.0*MeV};
    Gaudi::Property<double> m_logWeightThres{this, "logWeightThres", 4.2};
    DataHandle<eic::RawCalorimeterHitCollection>
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
	    const auto &rawhits = *m_inputHitCollection.get();
        // Create output collections
        auto &clusters = *m_outputClusterCollection.createAndPut();

        info() << "we have " << rawhits.size() << " raw hits" << endmsg;

        // energy time reconstruction
        eic::CalorimeterHitCollection hits;
        for (auto &rh : rawhits) {
            float energy = rh.amplitude()/100.*MeV;
            if (energy >= m_minModuleEdep) {
                float time = rh.timeStamp();
                auto pos = m_geoSvc->cellIDPositionConverter()->position(rh.cellID0());
                hits.push_back(eic::CalorimeterHit{
                    rh.cellID0(), rh.cellID1(), energy, time, {pos.X(), pos.Y(), pos.Z()}, 0
                });
            }
        }
        info() << "we have " << hits.size() << " hits" << endmsg;

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

        // reconstruct hit position for the cluster
        for (auto &cl : clusters) {
            reconstruct(cl);
            info() << cl.energy()/GeV << " GeV, (" << cl.position()[0]/mm << ", "
                   << cl.position()[1]/mm << ", " << cl.position()[2]/mm << ")" << endmsg;
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

        return (std::abs(pos1.x - pos2.x) <= (dim1[0] + dim2[0])/1.3) &&
               (std::abs(pos1.y - pos2.y) <= (dim1[1] + dim2[1])/1.3);
    }

    // grouping function with Depth-First Search
    void dfs_group(eic::Cluster group, int idx, eic::CalorimeterHitCollection &hits, std::vector<bool> &visits)
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
            auto cl = clusters.create();
            for (auto hit : group.hits()) {
                auto shit = hit.clone();
                cl.addhits(shit);
                scoll.push_back(shit);
            }
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

    void reconstruct(eic::Cluster cl) {
        float totalE = 0.;
        for (auto &hit : cl.hits()) {
            totalE += hit.energy();
        }
        cl.energy(totalE);

        // center of gravity with logarithmic weighting
        float totalW = 0., x = 0., y = 0., z = 0.;
        for (auto &hit : cl.hits()) {
            float weight = m_logWeightThres + std::log(hit.energy()/totalE);
            totalW += weight;
            x += hit.position().x * weight;
            y += hit.position().y * weight;
            z += hit.position().z * weight;
        }
        cl.position() = std::array<float, 3>{x/totalW, y/totalW, z/totalW};
    }
};

DECLARE_COMPONENT(CalorimeterIslandCluster)

} // namespace Jug::Reco

