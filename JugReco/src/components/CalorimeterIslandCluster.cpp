#include <algorithm>

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Property.h"
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
/*  Island Clustering Algorithm for Calorimeter Blocks
 *  1. group all the adjacent modules with the energy deposit above <minModuleEdep>
 *  2. split the groups between their local maxima with the energy deposit above <minClusterCenterEdep>
 *  3. reconstruct the clustrers
 */
class CalorimeterIslandCluster : public GaudiAlgorithm
{
public:
    Gaudi::Property<double> m_minModuleEdep{this, "minModuleEdep", 5.0*MeV};
    Gaudi::Property<double> m_minClusterCenterEdep{this, "minClusterCenterEdep", 50.0*MeV};
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
        auto clusterhits = m_outputClusterCollection.createAndPut();

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

        // group neighboring hits
        std::vector<bool> visits(hits.size(), false);
        std::vector<std::vector<eic::CalorimeterHit>> groups;
        for(size_t i = 0; i < hits.size(); ++i)
        {
            // already in a group
            if (visits[i]) {
                continue;
            }

            // create a new group and reserve some space for the possible hits
            groups.emplace_back();
            // group all the possible hits
            dfs_group(groups.back(), i, hits, visits);
        }
        return StatusCode::SUCCESS;
    }

private:
    // helper function to group hits
    inline bool is_neighbor(const eic::CalorimeterHit &h1, const eic::CalorimeterHit &h2)
    {
        auto pos1 = h1.position();
        auto pos2 = h2.position();
        auto dim1 = m_geoSvc->cellIDPositionConverter()->cellDimensions(h1.cellID0());
        auto dim2 = m_geoSvc->cellIDPositionConverter()->cellDimensions(h2.cellID0());

        return (std::abs(pos1.x - pos2.x) <= (dim1[0] + dim2[0])/2.) &&
               (std::abs(pos1.y - pos2.y) <= (dim1[1] + dim2[1])/2.);
    }

    // recursive function for the DFS grouping
    void dfs_group(std::vector<eic::CalorimeterHit> &group, int idx,
                   eic::CalorimeterHitCollection &hits, std::vector<bool> &visits)
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

DECLARE_COMPONENT(CalorimeterIslandCluster)

} // namespace Jug::Reco

