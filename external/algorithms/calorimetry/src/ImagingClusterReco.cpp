// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Chao Peng, Wouter Deconinck

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
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/MCRecoClusterParticleAssociationCollection.h"
#include "edm4eic/ProtoClusterCollection.h"
#include "edm4eic/vector_utils.h"

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
private:
  Gaudi::Property<int> m_trackStopLayer{this, "trackStopLayer", 9};

  DataHandle<edm4eic::ProtoClusterCollection> m_inputProtoClusters{"inputProtoClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::ClusterCollection> m_outputLayers{"outputLayers", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4eic::ClusterCollection> m_outputClusters{"outputClusters", Gaudi::DataHandle::Reader, this};

  // Collection for MC hits when running on MC
  Gaudi::Property<std::string> m_mcHits{this, "mcHits", ""};
  // Optional handle to MC hits
  std::unique_ptr<DataHandle<edm4hep::SimCalorimeterHitCollection>> m_mcHits_ptr;

  // Collection for associations when running on MC
  Gaudi::Property<std::string> m_outputAssociations{this, "outputAssociations", ""};
  // Optional handle to MC hits
  std::unique_ptr<DataHandle<edm4eic::MCRecoClusterParticleAssociationCollection>> m_outputAssociations_ptr;

public:
  ImagingClusterReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputProtoClusters", m_inputProtoClusters, "");
    declareProperty("outputLayers", m_outputLayers, "");
    declareProperty("outputClusters", m_outputClusters, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    // Initialize the optional MC input hit collection if requested
    if (m_mcHits != "") {
      m_mcHits_ptr =
        std::make_unique<DataHandle<edm4hep::SimCalorimeterHitCollection>>(m_mcHits, Gaudi::DataHandle::Reader,
        this);
    }

    // Initialize the optional association collection if requested
    if (m_outputAssociations != "") {
      m_outputAssociations_ptr =
        std::make_unique<DataHandle<edm4eic::MCRecoClusterParticleAssociationCollection>>(m_outputAssociations, Gaudi::DataHandle::Writer,
        this);
    }

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& proto = *m_inputProtoClusters.get();
    // output collections
    auto& layers   = *m_outputLayers.createAndPut();
    auto& clusters = *m_outputClusters.createAndPut();

    // Optional input MC data
    const edm4hep::SimCalorimeterHitCollection* mcHits = nullptr;
    if (m_mcHits_ptr) {
      mcHits = m_mcHits_ptr->get();
    }

    // Optional output associations
    edm4eic::MCRecoClusterParticleAssociationCollection* associations = nullptr;
    if (m_outputAssociations_ptr) {
      associations = m_outputAssociations_ptr->createAndPut();
    }

    for (const auto& pcl : proto) {
      if (!pcl.getHits().empty() && !pcl.getHits(0).isAvailable()) {
        warning() << "Protocluster hit relation is invalid, skipping protocluster" << endmsg;
        continue;
      }
      // get cluster and associated layers
      auto cl        = reconstruct_cluster(pcl);
      auto cl_layers = reconstruct_cluster_layers(pcl);

      // Get cluster direction from the layer profile
      auto [theta, phi] = fit_track(cl_layers);
      cl.setIntrinsicTheta(theta);
      cl.setIntrinsicPhi(phi);
      // no error on the intrinsic direction TODO

      // store layer and clusters on the datastore
      for (auto& layer : cl_layers) {
        layers.push_back(layer);
        cl.addToClusters(layer);
      }
      clusters.push_back(cl);


      // If mcHits are available, associate cluster with MCParticle
      if (m_mcHits_ptr.get() != nullptr && m_outputAssociations_ptr.get() != nullptr) {

        // 1. find pclhit with largest energy deposition
        auto pclhits = pcl.getHits();
        auto pclhit = std::max_element(
          pclhits.begin(),
          pclhits.end(),
          [](const auto& pclhit1, const auto& pclhit2) {
            return pclhit1.getEnergy() < pclhit2.getEnergy();
          }
        );

        // 2. find mchit with same CellID
        auto mchit = mcHits->begin();
        for ( ; mchit != mcHits->end(); ++mchit) {
          // break loop when CellID match found
          if (mchit->getCellID() == pclhit->getCellID()) {
            break;
          }
        }
        if (!(mchit != mcHits->end())) {
          // break if no matching hit found for this CellID
          warning() << "Proto-cluster has highest energy in CellID " << pclhit->getCellID()
                    << ", but no mc hit with that CellID was found." << endmsg;
          break;
        }

        // 3. find mchit's MCParticle
        const auto& mcp = mchit->getContributions(0).getParticle();

        // set association
        edm4eic::MutableMCRecoClusterParticleAssociation clusterassoc;
        clusterassoc.setRecID(cl.getObjectID().index);
        clusterassoc.setSimID(mcp.getObjectID().index);
        clusterassoc.setWeight(1.0);
        clusterassoc.setRec(cl);
        //clusterassoc.setSim(mcp);
        associations->push_back(clusterassoc);
      }

    }

    // debug output
    if (msgLevel(MSG::DEBUG)) {
      for (const auto& cl : clusters) {
        debug() << fmt::format("Cluster {:d}: Edep = {:.3f} MeV, Dir = ({:.3f}, {:.3f}) deg", cl.id(),
                               cl.getEnergy() * 1000., cl.getIntrinsicTheta() / M_PI * 180.,
                               cl.getIntrinsicPhi() / M_PI * 180.)
                << endmsg;
      }
    }

    return StatusCode::SUCCESS;
  }

private:
  template <typename T> static inline T pow2(const T& x) { return x * x; }

  static std::vector<edm4eic::Cluster> reconstruct_cluster_layers(const edm4eic::ProtoCluster& pcl) {
    const auto& hits    = pcl.getHits();
    const auto& weights = pcl.getWeights();
    // using map to have hits sorted by layer
    std::map<int, std::vector<std::pair<edm4eic::CalorimeterHit, float>>> layer_map;
    for (unsigned i = 0; i < hits.size(); ++i) {
      const auto hit = hits[i];
      auto lid       = hit.getLayer();
      if (layer_map.count(lid) == 0) {
        layer_map[lid] = {};
      }
      layer_map[lid].push_back({hit, weights[i]});
    }

    // create layers
    std::vector<edm4eic::Cluster> cl_layers;
    for (const auto& [lid, layer_hits] : layer_map) {
      auto layer = reconstruct_layer(layer_hits);
      cl_layers.push_back(layer);
    }
    return cl_layers;
  }

  static edm4eic::Cluster reconstruct_layer(const std::vector<std::pair<edm4eic::CalorimeterHit, float>>& hits) {
    edm4eic::MutableCluster layer;
    layer.setType(ClusterType::kClusterSlice);
    // Calculate averages
    double energy{0};
    double energyError{0};
    double time{0};
    double timeError{0};
    double sumOfWeights{0};
    auto pos            = layer.getPosition();
    for (const auto& [hit, weight] : hits) {
      energy += hit.getEnergy() * weight;
      energyError += std::pow(hit.getEnergyError() * weight, 2);
      time += hit.getTime() * weight;
      timeError += std::pow(hit.getTimeError() * weight, 2);
      pos = pos + hit.getPosition() * weight;
      sumOfWeights += weight;
      layer.addToHits(hit);
    }
    layer.setEnergy(energy);
    layer.setEnergyError(std::sqrt(energyError));
    layer.setTime(time / sumOfWeights);
    layer.setTimeError(std::sqrt(timeError) / sumOfWeights);
    layer.setNhits(hits.size());
    layer.setPosition(pos / sumOfWeights);
    // positionError not set
    // Intrinsic direction meaningless in a cluster layer --> not set

    // Calculate radius as the standard deviation of the hits versus the cluster center
    double radius = 0.;
    for (const auto& [hit, weight] : hits) {
      radius += std::pow(edm4eic::magnitude(hit.getPosition() - layer.getPosition()), 2);
    }
    layer.addToShapeParameters(std::sqrt(radius / layer.getNhits()));
    // TODO Skewedness

    return layer;
  }

  edm4eic::MutableCluster reconstruct_cluster(const edm4eic::ProtoCluster& pcl) {
    edm4eic::MutableCluster cluster;

    const auto& hits    = pcl.getHits();
    const auto& weights = pcl.getWeights();

    cluster.setType(ClusterType::kCluster3D);
    double energy      = 0.;
    double energyError = 0.;
    double time        = 0.;
    double timeError   = 0.;
    double meta        = 0.;
    double mphi        = 0.;
    double r           = 9999 * cm;
    for (unsigned i = 0; i < hits.size(); ++i) {
      const auto& hit    = hits[i];
      const auto& weight = weights[i];
      energy += hit.getEnergy() * weight;
      energyError += std::pow(hit.getEnergyError() * weight, 2);
      // energy weighting for the other variables
      const double energyWeight = hit.getEnergy() * weight;
      time += hit.getTime() * energyWeight;
      timeError += std::pow(hit.getTimeError() * energyWeight, 2);
      meta += edm4eic::eta(hit.getPosition()) * energyWeight;
      mphi += edm4eic::angleAzimuthal(hit.getPosition()) * energyWeight;
      r = std::min(edm4eic::magnitude(hit.getPosition()), r);
      cluster.addToHits(hit);
    }
    cluster.setEnergy(energy);
    cluster.setEnergyError(std::sqrt(energyError));
    cluster.setTime(time / energy);
    cluster.setTimeError(std::sqrt(timeError) / energy);
    cluster.setNhits(hits.size());
    cluster.setPosition(edm4eic::sphericalToVector(r, edm4eic::etaToAngle(meta / energy), mphi / energy));

    // shower radius estimate (eta-phi plane)
    double radius = 0.;
    for (const auto& hit : hits) {
      radius += pow2(edm4eic::eta(hit.getPosition()) - edm4eic::eta(cluster.getPosition())) +
                pow2(edm4eic::angleAzimuthal(hit.getPosition()) - edm4eic::angleAzimuthal(cluster.getPosition()));
    }
    cluster.addToShapeParameters(std::sqrt(radius / cluster.getNhits()));
    // Skewedness not calculated TODO

    // Optionally store the MC truth associated with the first hit in this cluster
    // FIXME no connection between cluster and truth in edm4hep
    // if (mcHits) {
    //  const auto& mc_hit    = (*mcHits)[pcl.getHits(0).ID.value];
    //  cluster.mcID({mc_hit.truth().trackID, m_kMonteCarloSource});
    //}

    return cluster;
  }

  std::pair<double /* polar */, double /* azimuthal */> fit_track(const std::vector<edm4eic::Cluster>& layers) const {
    int nrows = 0;
    decltype(edm4eic::ClusterData::position) mean_pos{0, 0, 0};
    for (const auto& layer : layers) {
      if ((layer.getNhits() > 0) && (layer.getHits(0).getLayer() <= m_trackStopLayer)) {
        mean_pos = mean_pos + layer.getPosition();
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
      if ((layer.getNhits() > 0) && (layer.getHits(0).getLayer() <= m_trackStopLayer)) {
        auto delta = layer.getPosition() - mean_pos;
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

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ImagingClusterReco)

} // namespace Jug::Reco
