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
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "eicd/CalorimeterHitCollection.h"
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

  // Collection for MC hits when running on MC
  Gaudi::Property<std::string> m_mcHits{this, "mcHits", ""};
  // Monte Carlo particle source identifier
  const int32_t m_kMonteCarloSource{uniqueID<int32_t>("MCParticles")};
  // Optional handle to MC hits
  std::unique_ptr<DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputMC;

  ImagingClusterReco(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin<>(name, info()) {
    declareProperty("inputProtoClusterCollection", m_inputProtoClusterCollection, "");
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputLayerCollection", m_outputLayerCollection, "");
    declareProperty("outputClusterCollection", m_outputClusterCollection, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    // Initialize the MC input hit collection if requested
    if (m_mcHits != "") {
      m_inputMC =
          std::make_unique<DataHandle<edm4hep::SimCalorimeterHitCollection>>(m_mcHits, Gaudi::DataHandle::Reader, this);
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
    // Optional MC data
    const edm4hep::SimCalorimeterHitCollection* mcHits = nullptr;
    if (m_inputMC) {
      mcHits = m_inputMC->get();
    }

    for (const auto& pcl : proto) {
      // get cluster and associated layers
      auto cl        = reconstruct_cluster(pcl, hits, mcHits);
      auto cl_layers = reconstruct_cluster_layers(pcl, hits);

      // Get cluster direction from the layer profile
      cl.direction(fit_track(cl_layers, m_trackStopLayer));

      // store layer and clusters on the datastore
      for (auto& layer : cl_layers) {
        // unique ID for this clusterlayer, starting at the last
        // cluster ID value to guarantee uniqueness
        layer.ID({static_cast<int32_t>(proto.size() + layers.size()), algorithmID()});
        layers.push_back(layer);
      }
      clusters.push_back(cl);
    }

    // debug output
    if (msgLevel(MSG::DEBUG)) {
      for (const auto& cl : clusters) {
        debug() << fmt::format("Cluster {:d}: Edep = {:.3f} MeV, Dir = ({:.3f}, {:.3f}) deg", cl.ID().value,
                               cl.energy() * 1000., cl.direction().theta / M_PI * 180., cl.direction().phi / M_PI * 180.)
                << endmsg;
      }
    }

    return StatusCode::SUCCESS;
  }

private:
  template <typename T> static inline T pow2(const T& x) { return x * x; }

  std::vector<eic::ClusterLayer> reconstruct_cluster_layers(const eic::ConstProtoCluster& pcl,
                                                            const eic::CalorimeterHitCollection& hits) const {
    // using map to have hits sorted by layer
    std::map<int, std::vector<std::pair<eic::ConstCalorimeterHit, eic::Weight>>> layer_map;
    for (const auto& clhit : pcl.hits()) {
      const auto hit = hits[clhit.index];
      auto lid = hit.layer();
      if (!layer_map.count(lid)) {
        layer_map[lid] = {};
      }
      layer_map[lid].push_back({hit, clhit.weight});
    }

    // create layers
    std::vector<eic::ClusterLayer> cl_layers;
    for (const auto& [lid, layer_hits] : layer_map) {
      auto layer = reconstruct_layer(pcl.ID().value, lid, layer_hits);
      cl_layers.push_back(layer);
    }
    return cl_layers;
  }

  eic::ClusterLayer reconstruct_layer(const int cid, const int lid,
                                      const std::vector<std::pair<eic::ConstCalorimeterHit, eic::Weight>>& hits) const {
    // use full members initialization here so it could catch changes in eicd
    eic::ClusterLayer layer{{}, {cid, algorithmID()}, lid, static_cast<uint32_t>(hits.size()), 0., 0., 0., 0., {}};

    // mean position and total energy
    eic::VectorXYZ pos;
    double energy = 0.;
    for (const auto& [hit, weight] : hits) {
      pos   = pos.add(hit.position());
      energy += hit.energy() * weight;
    }

    pos = pos.scale(1 / layer.nhits());
    layer.position(pos);
    layer.energy(energy);

    double radius = 0.;
    for (const auto& [hit, weight] : hits) {
      radius += hit.position().subtract(layer.position()).mag();
    }
    layer.radius(radius / layer.nhits());
    return layer;
  }

  eic::Cluster reconstruct_cluster(const eic::ConstProtoCluster& pcl, const eic::CalorimeterHitCollection& hits,
                                   const edm4hep::SimCalorimeterHitCollection* mcHits) {
    eic::Cluster cluster;
    cluster.ID({pcl.ID(), algorithmID()});
    // eta, phi center, weighted by energy
    double meta = 0.;
    double mphi = 0.;
    double energy = 0.;
    float r     = 9999 * cm;
    for (const auto& clhit : pcl.hits()) {
      const auto& hit = hits[clhit.index];
      meta += hit.position().eta() * hit.energy() * clhit.weight;
      mphi += hit.position().phi() * hit.energy() * clhit.weight;
      energy += hit.energy() * clhit.weight;
      r = std::min(hit.position().r(), r);
    }
    const double eta   = meta / energy;
    const double phi   = mphi / energy;
    const double theta = 2. * std::atan(std::exp(-eta));
    cluster.nhits(pcl.hits_size());
    cluster.energy(energy / m_sampFrac); // simple energy reconstruction //DEPRECATED
    eic::VectorPolar polar{r, theta, phi > M_PI ? phi - M_PI : phi};
    cluster.position(polar);

    // shower radius estimate (eta-phi plane)
    double radius = 0.;
    for (const auto& clhit : pcl.hits()) {
      const auto& hit = hits[clhit.index];
      radius += std::sqrt(pow2(hit.position().eta() - cluster.position().eta()) +
                          pow2(hit.position().phi() - cluster.position().phi()));
    }
    cluster.radius(radius / cluster.nhits());

    // Optionally store the MC truth associated with the first hit in this cluster
    if (mcHits) {
      const auto& mc_hit    = (*mcHits)[pcl.hits(0).ID.value];
      cluster.mcID({mc_hit.truth().trackID, m_kMonteCarloSource});
    }
    
    // fill additional info fields;
    cluster.polar(cluster.position());
    cluster.eta(cluster.position().eta());

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

