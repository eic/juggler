/*
 *  Reconstruct the cluster for imaging calorimeter
 *  Logarithmic weighting is used for mimicing energy deposit in transverse direction
 *
 *  Author: Chao Peng (ANL), 09/27/2020
 */
#include <algorithm>
#include <Eigen/Dense>
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

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/Utilities/Utils.hpp"

// Event Model related classes
#include "eicd/CalorimeterHit.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ImagingLayerCollection.h"
#include "eicd/ImagingClusterCollection.h"

using namespace Gaudi::Units;
using namespace Eigen;

namespace Jug::Reco {
class ImagingClusterReco : public GaudiAlgorithm
{
public:
    Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.6};
    Gaudi::Property<double> m_minEdep{this, "mininumEnergyDeposit", 0.5*MeV};
    Gaudi::Property<int> m_minNhits{this, "miniumNumberOfHits", 10};
    Gaudi::Property<int> m_trackStopLayer{this, "trackStopLayer", 9};
    Gaudi::Property<std::vector<int>> u_layerIDMaskRange{this, "layerIDMaskRange", {}};

    DataHandle<eic::ClusterCollection>
        m_inputClusterCollection{"inputClusterCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ImagingClusterCollection>
        m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer, this};
    DataHandle<eic::ImagingLayerCollection>
        m_outputLayerCollection{"outputLayerCollection", Gaudi::DataHandle::Writer, this};
    // Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;
    double m_depthCorr;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ImagingClusterReco(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputClusterCollection", m_inputClusterCollection, "");
        declareProperty("outputClusterCollection", m_outputClusterCollection, "");
        declareProperty("outputLayerCollection", m_outputLayerCollection, "");
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
        auto &input = *m_inputClusterCollection.get();
        auto &output = *m_outputClusterCollection.createAndPut();
        auto &layers = *m_outputLayerCollection.createAndPut();

        int ncl = 0;
        for (auto cl : input) {
            // basic information inherited from the cluster
            eic::ImagingCluster img;
            double edep = 0.;
            for (auto hit : cl.hits()) {
                edep += hit.energy();
            }
            img.nhits(cl.hits().size());
            img.edep(edep);
            img.energy(edep / m_sampFrac);

            // check if the cluster passes the filter
            if ((img.nhits() < m_minNhits) || (img.edep() < m_minEdep)) {
                continue;
            }

            // group hits to layers
            group_by_layer(img, cl, layers, ncl++);

            // fit intrinsic theta/phi
            fit_track(img, m_trackStopLayer);
            output.push_back(img);
        }

        for (auto [k, cl] : Jug::Utils::Enumerate(output)) {
            debug() << fmt::format("Cluster {:d}: Edep = {:.3f} MeV, Dir = ({:.3f}, {:.3f}) deg",
                                   k + 1, cl.edep()/MeV, cl.cl_theta()/M_PI*180., cl.cl_phi()/M_PI*180.)
                    << endmsg;
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

    void group_by_layer(eic::ImagingCluster &img, eic::Cluster &cluster, eic::ImagingLayerCollection &container, int cid)
    const
    {
        // using map to have id sorted
        std::map<int, std::vector<int>> hits_map;

        // group hits
        for (auto [ih, hit] : Jug::Utils::Enumerate(cluster.hits())) {
            auto lid = get_subid(hit.cellID(), id_mask, id_shift);
            auto it = hits_map.find(lid);
            if (it == hits_map.end()) {
                hits_map[lid] = {ih};
            } else {
                it->second.push_back(ih);
            }
        }

        // create layers
        for (auto it : hits_map) {
            eic::ImagingLayer layer;
            layer.clusterID(cid);
            layer.layerID(it.first);
            layer.edep(0.);
            layer.position({0., 0., 0.});
            double mx = 0., my = 0., mz = 0.;
            for (auto [k, hid] : Jug::Utils::Enumerate(it.second)) {
                if (k >= layer.hits().size()) {
                    warning() << fmt::format("Discard hit {:d} at layer {:d} because container caps at {:d}",
                                             k + 1, it.first, layer.hits().size())
                              << endmsg;
                    continue;
                }
                auto hit = cluster.hits(hid);
                layer.hits()[k] = eic::BaseHit{hit.x(), hit.y(), hit.z(), hit.energy()};
                mx += hit.x();
                my += hit.y();
                mz += hit.z();
                layer.edep(layer.edep() + hit.energy());
                layer.nhits(k + 1);
            }
            layer.x(mx/layer.nhits());
            layer.y(my/layer.nhits());
            layer.z(mz/layer.nhits());
            // add relation
            container.push_back(layer);
            img.addlayers(layer);
        }
    }

    void fit_track(eic::ImagingCluster &img, int stop_layer) const
    {
        int nrows = 0;
        double mx = 0., my = 0., mz = 0.;
        for (auto layer : img.layers()) {
            if ((layer.layerID() <= stop_layer) && (layer.nhits() > 0)) {
                mx += layer.x();
                my += layer.y();
                mz += layer.z();
                nrows ++;
            }
        }
        // cannot fit
        if (nrows < 2) {
            return;
        }

        mx /= nrows;
        my /= nrows;
        mz /= nrows;
        // fill position data
        MatrixXd pos(nrows, 3);
        int ir = 0;
        for (auto layer : img.layers()) {
            if ((layer.layerID() <= stop_layer) && (layer.nhits() > 0)) {
                pos(ir, 0) = layer.x() - mx;
                pos(ir, 1) = layer.y() - my;
                pos(ir, 2) = layer.z() - mz;
                ir ++;
            }
        }

        JacobiSVD<MatrixXd> svd(pos, ComputeThinU | ComputeThinV);
        // debug() << pos << endmsg;
        // debug() << svd.matrixV() << endmsg;
        auto dir = svd.matrixV().col(0);
        img.cl_theta(std::acos(dir(2)));
        img.cl_phi(std::atan2(dir(1), dir(0)));
        // extract 3d line with SVD
    }

};

DECLARE_COMPONENT(ImagingClusterReco)

} // namespace Jug::Reco

