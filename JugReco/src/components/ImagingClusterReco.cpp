/*
 *  Reconstruct the cluster/layer info for imaging calorimeter
 *  Logarithmic weighting is used to describe energy deposit in transverse direction
 *
 *  Author: Chao Peng (ANL), 06/02/2021
 */
#include "fmt/format.h"
#include <Eigen/Dense>
#include <algorithm>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/UniqueID.h"
#include "JugBase/Utilities/Utils.hpp"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/Cluster3DInfoCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ClusterLayerCollection.h"
#include "eicd/ProtoClusterCollection.h"
#include "eicd/VectorPolar.h"
#include "eicd/VectorXYZ.h"

using namespace Gaudi::Units;
using namespace Eigen;

namespace Jug::Reco {

/** Imaging cluster reconstruction.
 *
 *  Reconstruct the cluster/layer info for imaging calorimeter
 *  Logarithmic weighting is used to describe energy deposit in transverse direction
 *
 *  \ingroup reco
 */
class ImagingClusterReco : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};
  Gaudi::Property<int> m_trackStopLayer{this, "trackStopLayer", 9};

  DataHandle<eic::ProtoClusterCollection> m_inputProtoClusterCollection{"inputProtoClusterCollection",
                                                                        Gaudi::DataHandle::Reader, this};
  DataHandle<eic::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ClusterLayerCollection> m_outputLayerCollection{"outputLayerCollection", Gaudi::DataHandle::Writer,
                                                                  this};
  DataHandle<eic::ClusterCollection> m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Reader,
                                                               this};
  DataHandle<eic::Cluster3DInfoCollection> m_outputInfoCollection{"outputInfoCollection", Gaudi::DataHandle::Reader,
                                                                  this};

  ImagingClusterReco(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin<>(name, info()) {
    declareProperty("inputProtoClusterCollection", m_inputProtoClusterCollection, "");
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputLayerCollection", m_outputLayerCollection, "");
    declareProperty("outputClusterCollection", m_outputClusterCollection, "");
    declareProperty("outputInfoCollection", m_outputInfoCollection, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& proto = *m_inputProtoClusterCollection.get();
    const auto& hits  = *m_inputHitCollection.get();
    // output collections
    auto& layers   = *m_outputLayerCollection.createAndPut();
    auto& clusters = *m_outputClusterCollection.createAndPut();
    auto& info     = *m_outputInfoCollection.createAndPut();

    // Create a map of clusterID --> associated ProtoCluster by looping over our clustered hits
    std::map<int, std::vector<std::pair<eic::ConstProtoCluster, 
                                        eic::ConstCalorimeterHit>>> cluster_map;
    for (const auto& pc : proto) {
      const size_t clusterID = pc.clusterID();
      if (!cluster_map.count(clusterID)) {
        cluster_map[clusterID] = {};
      }
      size_t idx;
      for (idx = 0; idx < hits.size(); ++idx) {
        if (hits[idx].ID() == pc.hitID()) {
          break;
        }
      }
      cluster_map[clusterID].push_back({pc, hits[idx]});
    }

    for (const auto& [cid, hit_info] : cluster_map) {
      // get cluster and associated layers
      auto cl        = reconstruct_cluster(hit_info, cid);
      auto cl_layers = reconstruct_cluster_layers(hit_info, cid);

      // reconstruct cluster direction
      eic::Cluster3DInfo cl_info{cl.ID(), cl.position(), cl.position().eta(), fit_track(cl_layers, m_trackStopLayer)};

      // store layer and clusters on the datastore
      for (auto& layer : cl_layers) {
        layer.ID(layers.size()); // unique ID for this clusterlayer
        layers.push_back(layer);
        // cl.addlayers(layer); // deprectated
      }
      clusters.push_back(cl);
      info.push_back(cl_info);
    }

    // debug output
    int idx = 0;
    if (msgLevel(MSG::DEBUG)) {
      for (const auto& cl : clusters) {
        debug() << fmt::format("Cluster {:d}: Edep = {:.3f} MeV, Dir = ({:.3f}, {:.3f}) deg", cl.ID(),
                               cl.energy() * 1000., info[idx].direction().theta / M_PI * 180.,
                               info[idx].direction().phi / M_PI * 180.)
                << endmsg;
        idx += 1;
      }
    }

    return StatusCode::SUCCESS;
  }

private:
  template <typename T> static inline T pow2(const T& x) { return x * x; }

  std::vector<eic::ClusterLayer>
  reconstruct_cluster_layers(const std::vector<std::pair<eic::ConstProtoCluster, eic::ConstCalorimeterHit>>& hit_info,
                             const int cid) const {
    // using map to have hits sorted by layer
    std::map<int, std::vector<std::pair<eic::ConstProtoCluster, 
                                        eic::ConstCalorimeterHit>>> layer_map;
    for (const auto& [proto, hit] : hit_info) {
      auto lid = hit.layer();
      if (!layer_map.count(lid)) {
        layer_map[lid] = {};
      }
      layer_map[lid].push_back({proto, hit});
    }

    // create layers
    std::vector<eic::ClusterLayer> cl_layers;
    for (const auto& [lid, layer_hit_info] : layer_map) {
      auto layer = reconstruct_layer(layer_hit_info, cid, lid);
      cl_layers.push_back(layer);
    }
    return cl_layers;
  }

  eic::ClusterLayer reconstruct_layer(const std::vector<std::pair<eic::ConstProtoCluster, eic::ConstCalorimeterHit>>& hit_info,
                                      const int cid, const int lid) const {
    // use full members initialization here so it could catch changes in ecid
    eic::ClusterLayer layer{-1, cid, lid, static_cast<uint32_t>(hit_info.size()), algorithmID(), 0., 0., 0., 0., {}};

    // mean position and total energy
    eic::VectorXYZ pos;
    double energy = 0.;
    for (const auto& [proto, hit] : hit_info) {
      pos   = pos.add(hit.position());
      energy += hit.energy() * proto.weight();
    }

    pos = pos.scale(1 / layer.nhits());
    layer.energy(energy);

    double radius = 0.;
    for (const auto& [proto, hit] : hit_info) {
      radius += std::sqrt(pow2(hit.position().x - layer.position().x) + pow2(hit.position().y - layer.position().y) +
                          pow2(hit.position().z - layer.position().z));
    }
    layer.radius(radius / layer.nhits());
    return layer;
  }

  eic::Cluster
  reconstruct_cluster(const std::vector<std::pair<eic::ConstProtoCluster, eic::ConstCalorimeterHit>>& hit_info,
                      const int cid) const {
    eic::Cluster cluster;
    cluster.ID(cid);
    cluster.source(algorithmID());
    // eta, phi center, weighted by energy
    double meta = 0.;
    double mphi = 0.;
    double energy = 0.;
    float r     = 9999 * cm;
    for (const auto& [proto, hit] : hit_info) {
      meta += hit.position().eta() * hit.energy() * proto.weight();
      mphi += hit.position().phi() * hit.energy() * proto.weight();
      energy += hit.energy() * proto.weight();
      r = std::min(hit.position().r(), r);
    }
    const double eta   = meta / energy;
    const double phi   = mphi / energy;
    const double theta = 2. * std::atan(std::exp(-eta));
    cluster.nhits(hit_info.size());
    cluster.energy(energy / m_sampFrac); // simple energy reconstruction //DEPRECATED
    eic::VectorPolar polar{r, theta, phi > M_PI ? phi - M_PI : phi};
    cluster.position(polar);

    // shower radius estimate (eta-phi plane)
    double radius = 0.;
    for (const auto& [proto, hit] : hit_info) {
      radius += std::sqrt(pow2(hit.position().eta() - cluster.position().eta()) +
                          pow2(hit.position().phi() - cluster.position().phi()));
    }
    cluster.radius(radius / cluster.nhits());

    return cluster;
  }

  eic::Direction fit_track(const std::vector<eic::ClusterLayer>& layers, const int stop_layer) const {
    int nrows = 0;
    double mx = 0.;
    double my = 0.;
    double mz = 0.;
    for (const auto& layer : layers) {
      if ((layer.layer() <= stop_layer) && (layer.nhits() > 0)) {
        mx += layer.position().x;
        my += layer.position().y;
        mz += layer.position().z;
        nrows += 1;
      }
    }
    // cannot fit
    if (nrows < 2) {
      return {};
    }

    mx /= nrows;
    my /= nrows;
    mz /= nrows;
    // fill position data
    MatrixXd pos(nrows, 3);
    int ir = 0;
    for (const auto& layer : layers) {
      if ((layer.layer() <= stop_layer) && (layer.nhits() > 0)) {
        pos(ir, 0) = layer.position().x - mx;
        pos(ir, 1) = layer.position().y - my;
        pos(ir, 2) = layer.position().z - mz;
        ir += 1;
      }
    }

    JacobiSVD<MatrixXd> svd(pos, ComputeThinU | ComputeThinV);
    // debug() << pos << endmsg;
    // debug() << svd.matrixV() << endmsg;
    const auto dir = svd.matrixV().col(0);
    // theta and phi
    return {std::acos(dir(2)), std::atan2(dir(1), dir(0))};
  }
};

DECLARE_COMPONENT(ImagingClusterReco)

} // namespace Jug::Reco

