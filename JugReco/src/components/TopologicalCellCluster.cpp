/*
 *  Topological Cell Clustering Algorithm for Sampling Calorimetry
 *  1. group all the adjacent modules
 *
 *  Author: Chao Peng (ANL), 04/06/2021
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
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  class TopologicalCellCluster : public GaudiAlgorithm {
  public:
    Gaudi::Property<int> m_adjLayerDiff{this, "adjLayerDiff", 1};
    Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
    Gaudi::Property<std::string> m_readout{this, "readoutClass", "EcalBarrelHits"};
    Gaudi::Property<std::string> m_layerField{this, "layerField", "layer"};
    Gaudi::Property<std::string> m_sectorField{this, "sectorField", "sector"};
    Gaudi::Property<std::vector<double>> u_localRanges{this, "localRanges", {1.0*mm, 1.0*mm}};
    Gaudi::Property<std::vector<double>> u_adjLayerRanges{this, "adjLayerRanges", {0.01*M_PI, 0.01*M_PI}};
    Gaudi::Property<double> m_minClusterCenterEdep{this, "minClusterCenterEdep", 0.5*MeV};
    DataHandle<eic::CalorimeterHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ClusterCollection>
        m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer, this};
    DataHandle<eic::CalorimeterHitCollection>
        m_splitHitCollection{"splitHitCollection", Gaudi::DataHandle::Writer, this};
    SmartIF<IGeoSvc> m_geoSvc;
    // visit readout fields
    dd4hep::BitFieldCoder *id_dec;
    size_t sector_idx, layer_idx;

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

        // check geometry service
        m_geoSvc = service(m_geoSvcName);
        if (!m_geoSvc) {
            error() << "Unable to locate Geometry Service. "
                    << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
            return StatusCode::FAILURE;
        }

        try {
            id_dec = m_geoSvc->detector()->readout(m_readout).idSpec().decoder();
            sector_idx = id_dec->index(m_sectorField);
            layer_idx = id_dec->index(m_layerField);
        } catch (...) {
            error() << "Failed to load ID decoder for " << m_readout << endmsg;
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
            debug() << fmt::format("hit {}, layer {}, sector {}", i,
                                   id_dec->get(hits[i].cellID(), layer_idx),
                                   id_dec->get(hits[i].cellID(), sector_idx))
                    << endmsg;
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
    bool is_neighbor(const eic::ConstCalorimeterHit &h1, const eic::ConstCalorimeterHit &h2) const
    {
        // we will merge different sectors later using global positions
        int s1 = id_dec->get(h1.cellID(), sector_idx);
        int s2 = id_dec->get(h2.cellID(), sector_idx);
        if (s1 != s2) { return false; }

        int l1 = id_dec->get(h1.cellID(), layer_idx);
        int l2 = id_dec->get(h2.cellID(), layer_idx);

        // layer check
        int ldiff = std::abs(l1 - l2);
        // same layer, check local positions
        if (!ldiff) {
            return (std::abs(h1.local_x() - h2.local_x()) <= u_localRanges[0]) &&
                   (std::abs(h1.local_y() - h2.local_y()) <= u_localRanges[1]);
        } else if (ldiff <= m_adjLayerDiff) {
            double eta1, phi1, r1, eta2, phi2, r2;
            calc_eta_phi_r(h1.x(), h1.y(), h1.z(), eta1, phi1, r1);
            calc_eta_phi_r(h2.x(), h2.y(), h2.z(), eta2, phi2, r2);

            return (std::abs(eta1 - eta2) <= u_adjLayerRanges[0]) &&
                   (std::abs(phi1 - phi2) <= u_adjLayerRanges[1]);
        }

        // not in adjacent layers
        return false;
    }

    static void calc_eta_phi_r(double x, double y, double z, double &eta, double &phi, double &r) {
        r = std::sqrt(x*x + y*y + z*z);
        phi = std::atan2(y, x);
        double theta = std::acos(z/r);
        eta = -std::log(std::tan(theta/2.));
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

