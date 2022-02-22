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
#include "JugBase/Utilities/Utils.hpp"
#include "JugReco/ClusterTypes.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ProtoClusterCollection.h"
#include "eicd/vector_utils.h"

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
class ImagingClusterReco : public GaudiAlgorithm {
public:
  Gaudi::Property<int> m_trackStopLayer{this, "trackStopLayer", 9};

  DataHandle<eic::ProtoClusterCollection> m_inputProtoClusters{"inputProtoClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ClusterCollection> m_outputLayers{"outputLayers", Gaudi::DataHandle::Writer, this};
  DataHandle<eic::ClusterCollection> m_outputClusters{"outputClusters", Gaudi::DataHandle::Reader, this};

  // Collection for MC hits when running on MC
  Gaudi::Property<std::string> m_mcHits{this, "mcHits", ""};
  // Optional handle to MC hits
  std::unique_ptr<DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputMC;

  ImagingClusterReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputProtoClusters", m_inputProtoClusters, "");
    declareProperty("outputLayers", m_outputLayers, "");
    declareProperty("outputClusters", m_outputClusters, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    // TODO move to own algo
    // Initialize the MC input hit collection if requested
    // if (m_mcHits != "") {
    //  m_inputMC =
    //      std::make_unique<DataHandle<edm4hep::SimCalorimeterHitCollection>>(m_mcHits, Gaudi::DataHandle::Reader,
    //      this);
    //}

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& proto = *m_inputProtoClusters.get();
    // output collections
    auto& layers   = *m_outputLayers.createAndPut();
    auto& clusters = *m_outputClusters.createAndPut();
    // Optional MC data
    // TODO remove to other algo
    // const edm4hep::SimCalorimeterHitCollection* mcHits = nullptr;
    // if (m_inputMC) {
    //  mcHits = m_inputMC->get();
    //}

    for (const auto& pcl : proto) {
      if (!pcl.hits().empty() && !pcl.hits(0).isAvailable()) {
        warning() << "Protocluster hit relation is invalid, skipping protocluster" << endmsg;
        continue;
      }
      // get cluster and associated layers
      auto cl        = reconstruct_cluster(pcl);
      auto cl_layers = reconstruct_cluster_layers(pcl);

      // Get cluster direction from the layer profile
      auto [theta, phi] = fit_track(cl_layers);
      cl.intrinsicTheta(theta);
      cl.intrinsicPhi(phi);
      // no error on the intrinsic direction TODO

      // store layer and clusters on the datastore
      for (auto& layer : cl_layers) {
        layers.push_back(layer);
        cl.addclusters(layer);
      }
      clusters.push_back(cl);
    }

    // debug output
    if (msgLevel(MSG::DEBUG)) {
      for (const auto& cl : clusters) {
        debug() << fmt::format("Cluster {:d}: Edep = {:.3f} MeV, Dir = ({:.3f}, {:.3f}) deg", cl.id(),
                               cl.energy() * 1000., cl.intrinsicTheta() / M_PI * 180., cl.intrinsicPhi() / M_PI * 180.)
                << endmsg;
      }
    }

    return StatusCode::SUCCESS;
  }

private:
  template <typename T> static inline T pow2(const T& x) { return x * x; }

  std::vector<eic::Cluster> reconstruct_cluster_layers(const eic::ConstProtoCluster& pcl) const {
    const auto& hits    = pcl.hits();
    const auto& weights = pcl.weights();
    // using map to have hits sorted by layer
    std::map<int, std::vector<std::pair<eic::ConstCalorimeterHit, float>>> layer_map;
    for (unsigned i = 0; i < hits.size(); ++i) {
      const auto hit = hits[i];
      auto lid       = hit.layer();
      if (!layer_map.count(lid)) {
        layer_map[lid] = {};
      }
      layer_map[lid].push_back({hit, weights[i]});
    }

    // create layers
    std::vector<eic::Cluster> cl_layers;
    for (const auto& [lid, layer_hits] : layer_map) {
      auto layer = reconstruct_layer(layer_hits);
      cl_layers.push_back(layer);
    }
    return cl_layers;
  }

  eic::Cluster reconstruct_layer(const std::vector<std::pair<eic::ConstCalorimeterHit, float>>& hits) const {
    eic::Cluster layer;
    layer.type(ClusterType::kClusterSlice);
    // Calculate averages
    double energy;
    double energyError;
    double time;
    double timeError;
    double sumOfWeights = 0;
    auto pos            = layer.position();
    for (const auto& [hit, weight] : hits) {
      energy += hit.energy() * weight;
      energyError += std::pow(hit.energyError() * weight, 2);
      time += hit.time() * weight;
      timeError += std::pow(hit.timeError() * weight, 2);
      pos = pos + hit.position() * weight;
      sumOfWeights += weight;
      layer.addhits(hit);
    }
    layer.energy(energy);
    layer.energyError(std::sqrt(energyError));
    layer.time(time / sumOfWeights);
    layer.timeError(std::sqrt(timeError) / sumOfWeights);
    layer.nhits(hits.size());
    layer.position(pos / sumOfWeights);
    // positionError not set
    // Intrinsic direction meaningless in a cluster layer --> not set

    // Calculate radius as the standard deviation of the hits versus the cluster center
    double radius = 0.;
    for (const auto& [hit, weight] : hits) {
      radius += std::pow(eicd::magnitude(hit.position() - layer.position()), 2);
    }
    layer.addshapeParameters(std::sqrt(radius / layer.nhits()));
    // TODO Skewedness

    return layer;
  }

  eic::Cluster reconstruct_cluster(const eic::ConstProtoCluster& pcl) {
    eic::Cluster cluster;

    const auto& hits    = pcl.hits();
    const auto& weights = pcl.weights();

    cluster.type(ClusterType::kCluster3D);
    double energy      = 0.;
    double energyError = 0.;
    double time        = 0.;
    double timeError   = 0.;
    double meta        = 0.;
    double mphi        = 0.;
    double r           = 9999 * cm;
    for (unsigned i = 0; i < hits.size(); ++i) {
      const auto& hit   = hits[i];
      const auto weight = weights[i];
      energy += hit.energy() * weight;
      energyError += std::pow(hit.energyError() * weight, 2);
      // energy weighting for the other variables
      const double energyWeight = hit.energy() * weight;
      time += hit.time() * energyWeight;
      timeError += std::pow(hit.timeError() * energyWeight, 2);
      meta += eicd::eta(hit.position()) * energyWeight;
      mphi += eicd::angleAzimuthal(hit.position()) * energyWeight;
      r = std::min(eicd::magnitude(hit.position()), r);
      cluster.addhits(hit);
    }
    cluster.energy(energy);
    cluster.energyError(std::sqrt(energyError));
    cluster.time(time / energy);
    cluster.timeError(std::sqrt(timeError) / energy);
    cluster.nhits(hits.size());
    cluster.position(eicd::sphericalToVector(r, eicd::etaToAngle(meta / energy), mphi / energy));

    // shower radius estimate (eta-phi plane)
    double radius = 0.;
    for (const auto& hit : hits) {
      radius += pow2(eicd::eta(hit.position()) - eicd::eta(cluster.position())) +
                pow2(eicd::angleAzimuthal(hit.position()) - eicd::angleAzimuthal(cluster.position()));
    }
    cluster.addshapeParameters(std::sqrt(radius / cluster.nhits()));
    // Skewedness not calculated TODO

    // Optionally store the MC truth associated with the first hit in this cluster
    // FIXME no connection between cluster and truth in edm4hep
    // if (mcHits) {
    //  const auto& mc_hit    = (*mcHits)[pcl.hits(0).ID.value];
    //  cluster.mcID({mc_hit.truth().trackID, m_kMonteCarloSource});
    //}

    return cluster;
  }

  std::pair<double /* polar */, double /* azimuthal */> fit_track(const std::vector<eic::Cluster>& layers) const {
    int nrows = 0;
    eicd::Vector3f mean_pos{0, 0, 0};
    for (const auto& layer : layers) {
      if ((layer.nhits() > 0) && (layer.hits(0).layer() <= m_trackStopLayer)) {
        mean_pos = mean_pos + layer.position();
        nrows += 1;
      }
    }
    // cannot fit
    if (nrows < 2) {
      return {};
    }

    mean_pos = mean_pos / nrows;
    // fill position data
    MatrixXd pos(nrows, 3);
    int ir = 0;
    for (const auto& layer : layers) {
      if ((layer.nhits() > 0) && (layer.hits(0).layer() <= m_trackStopLayer)) {
        auto delta = layer.position() - mean_pos;
        pos(ir, 0) = delta.x;
        pos(ir, 1) = delta.y;
        pos(ir, 2) = delta.z;
        ir += 1;
      }
    }

    JacobiSVD<MatrixXd> svd(pos, ComputeThinU | ComputeThinV);
    const auto dir = svd.matrixV().col(0);
    // theta and phi
    return {std::acos(dir(2)), std::atan2(dir(1), dir(0))};
  }
};

DECLARE_COMPONENT(ImagingClusterReco)

} // namespace Jug::Reco

