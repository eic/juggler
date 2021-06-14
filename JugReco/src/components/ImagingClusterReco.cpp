/*
 *  Reconstruct the cluster/layer info for imaging calorimeter
 *  Logarithmic weighting is used to describe energy deposit in transverse direction
 *
 *  Author: Chao Peng (ANL), 06/02/2021
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
#include "eicd/ImagingPixel.h"
#include "eicd/ImagingLayerCollection.h"
#include "eicd/ImagingClusterCollection.h"

using namespace Gaudi::Units;
using namespace Eigen;

namespace Jug::Reco {

class ImagingClusterReco : public GaudiAlgorithm
{
public:
    Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.6};
    Gaudi::Property<int> m_trackStopLayer{this, "trackStopLayer", 9};

    DataHandle<eic::ImagingClusterCollection>
        m_inputClusterCollection{"inputClusterCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ImagingLayerCollection>
        m_outputLayerCollection{"outputLayerCollection", Gaudi::DataHandle::Writer, this};

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ImagingClusterReco(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputClusterCollection", m_inputClusterCollection, "");
        declareProperty("outputLayerCollection", m_outputLayerCollection, "");
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
        auto &input = *m_inputClusterCollection.get();
        auto &layers = *m_outputLayerCollection.createAndPut();

        int ncl = 1;
        for (auto cl : input) {
            // simple energy reconstruction
            cl.energy(cl.edep() / m_sampFrac);

            // group hits to layers
            group_by_layer(cl, layers, ncl++);

            for (auto &layer : layers) {
                reconstruct_layer(layer);
            }

            reconstruct_cluster(cl);
            // fit intrinsic theta/phi
            fit_track(cl, m_trackStopLayer);
        }

        for (auto [k, cl] : Jug::Utils::Enumerate(input)) {
            debug() << fmt::format("Cluster {:d}: Edep = {:.3f} MeV, Dir = ({:.3f}, {:.3f}) deg",
                                   k + 1, cl.edep()/MeV, cl.cl_theta()/M_PI*180., cl.cl_phi()/M_PI*180.)
                    << endmsg;
        }

        return StatusCode::SUCCESS;
    }

private:
    template<typename T> static inline T pow2(const T& x) { return x*x; }

    void group_by_layer(eic::ImagingCluster &cluster, eic::ImagingLayerCollection &container, int cid) const
    {
        // using map to have id sorted
        std::map<int, std::vector<size_t>> hits_map;

        // group hits
        for (auto [ih, hit] : Jug::Utils::Enumerate(cluster.hits())) {
            auto lid = hit.layerID();
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
            for (auto hid : it.second) {
                auto hit = cluster.hits(hid);
                layer.addhits(hit);
            }
            // add relation
            container.push_back(layer);
            cluster.addlayers(layer);
        }
    }

    void reconstruct_layer(eic::ImagingLayer layer) const
    {
        // number of hits
        layer.nhits(layer.hits_size());
        // mean position and total edep
        double mx = 0., my = 0., mz = 0., edep = 0.;
        for (auto hit : layer.hits()) {
            mx += hit.x();
            my += hit.y();
            mz += hit.z();
            edep += hit.edep();
        }

        layer.x(mx/layer.nhits());
        layer.y(my/layer.nhits());
        layer.z(mz/layer.nhits());
        layer.edep(edep);

        double radius = 0.;
        for (auto hit : layer.hits()) {
            radius += std::sqrt(pow2(hit.x() - layer.x()) + pow2(hit.y() - layer.y()) + pow2(hit.z() - layer.z()));
        }
        layer.radius(radius/layer.nhits());
    }

    void reconstruct_cluster(eic::ImagingCluster cluster) const
    {
        // eta, phi center, weighted by energy
        double meta = 0., mphi = 0., edep = 0.;
        for (auto hit : cluster.hits()) {
            meta += hit.eta()*hit.edep();
            mphi += hit.phi()*hit.edep();
            edep += hit.edep();
        }
        cluster.nhits(cluster.hits_size());
        cluster.edep(edep);
        cluster.eta(meta/edep);
        cluster.phi(mphi/edep);

        double radius = 0.;
        for (auto hit : cluster.hits()) {
            radius += std::sqrt(pow2(hit.eta() - cluster.eta()) + pow2(hit.phi() - cluster.phi()));
        }
        cluster.radius(radius/cluster.nhits());
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

