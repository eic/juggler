#include <limits>
#include <numbers>

#include <fmt/format.h>
// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"

#include "JugBase/DataHandle.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Fast {

/** Simple algorithm to merge the energy measurement from cluster1 with the position
 * measurement of cluster2 (in case matching clusters are found). If not, it will
 * propagate the raw cluster from cluster1 or cluster2
 *
 * Matching occurs based on the mc truth information of the clusters.
 *
 * \ingroup reco
 */
class TruthEnergyPositionClusterMerger : public GaudiAlgorithm {
public:
  // Input
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"MCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ClusterCollection> m_energyClusters{"EnergyClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ClusterCollection> m_positionClusters{"PositionClusters", Gaudi::DataHandle::Reader, this};
  // Output
  DataHandle<eic::ClusterCollection> m_outputClusters{"OutputClusters", Gaudi::DataHandle::Writer, this};

public:
  TruthEnergyPositionClusterMerger(const std::string& name, ISvcLocator* svcLoc)
<<<<<<< HEAD
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputMCParticles, "MCParticles");
=======
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticles, "mcparticles");
>>>>>>> Update for new EICD changes, with major changes to calorimetry
    declareProperty("inputEnergyClusters", m_energyClusters, "Cluster collection with good energy precision");
    declareProperty("inputPositionClusters", m_positionClusters, "Cluster collection with good position precision");
    declareProperty("outputClusters", m_outputClusters, "Cluster collection with good energy precision");
  }

  StatusCode initialize() override { return StatusCode::SUCCESS; }

  StatusCode execute() override {
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Merging energy and position clusters for new event" << endmsg;
    }
    // input
    const auto& mcparticles = *(m_inputMCParticles.get());
    const auto& energy_clus = *(m_energyClusters.get());
    const auto& pos_clus    = *(m_positionClusters.get());
    // output
    auto& merged    = *(m_outputClusters.createAndPut());

    if (!energy_clus.size() && !pos_clus.size()) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "Nothing to do for this event, returning..." << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 0/2: Getting indexed list of clusters..." << endmsg;
    }
    // get an indexed map of all clusters
    auto energyMap = indexedClusters(energy_clus);
    auto posMap    = indexedClusters(pos_clus);

    // loop over all position clusters and match with energy clusters
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 1/2: Matching all position clusters to the available energy clusters..." << endmsg;
    }
    for (const auto& [mcID, pclus] : posMap) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Processing position cluster " << pclus.id() << ", mcID: " << mcID << ", energy: " << pclus.energy()
                << endmsg;
      }
      if (energyMap.count(mcID)) {
        const auto& eclus = energyMap[mcID];
        auto new_clus = merged.create();
        new_clus.energy(eclus.energy());
        new_clus.energyError(eclus.energyError());
        new_clus.time(pclus.time());
        new_clus.nhits(pclus.nhits() + eclus.nhits());
        new_clus.position(pclus.position());
        new_clus.positionError(pclus.positionError());
        new_clus.mcID(mcID);
        new_clus.addclusters(pclus);
        new_clus.addclusters(eclus);
        for (const auto& cl : {pclus, eclus}) {
          for (const auto& hit : cl.hits()) {
            new_clus.addhits(hit);
            new_clus.addhitContributions(hit.energy());
          }
          new_clus.addsubdetectorEnergies(cl.energy());
        }
        for (const auto& param : pclus.shapeParameters()) {
          new_clus.addshapeParameters(param);
        }
        if (msgLevel(MSG::DEBUG)) {
          debug() << "   --> Found matching energy cluster " << eclus.id() << ", energy: " << eclus.energy() << endmsg;
          debug() << "   --> Created new combined cluster " << new_clus.id() << ", energy: " << new_clus.energy() << endmsg;
        }
        // erase the energy cluster from the map, so we can in the end account for all
        // remaining clusters
        energyMap.erase(mcID);
      } else {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "   --> No matching energy cluster found, copying over position cluster" << endmsg;
        }
        merged.push_back(pclus.clone());
        merged[merged.size()-1].addclusters(pclus);
      }
    }
    // Collect remaining energy clusters. Use mc truth position for these clusters, as
    // they should really have a match in the position clusters (and if they don't it due
    // to a clustering error).
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 2/2: Collecting remaining energy clusters..." << endmsg;
    }
    for (const auto& [mcID, eclus] : energyMap) {
      const auto& mc = mcparticles[mcID.value];
      const auto& p = mc.getMomentum();
      const auto phi = std::atan2(p.y, p.x);
      const auto theta = std::atan2(std::hypot(p.x, p.y), p.z);
      auto new_clus  = merged.create();
      new_clus.energy(eclus.energy());
      new_clus.energyError(eclus.energyError());
      new_clus.time(eclus.time());
      new_clus.nhits(eclus.nhits());
      // use nominal radius of 110cm, and use start vertex theta and phi
      new_clus.position(eic::VectorXYZ::fromSpherical(110, mc.ps().theta(), mc.ps().phi()));
      new_clus.mcID(mcID);
      new_clus.addclusters(eclus);
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Processing energy cluster " << eclus.id() << ", mcID: " << mcID << ", energy: " << eclus.energy()
                << endmsg;
        debug() << "   --> Created new 'combined' cluster " << new_clus.id() << ", energy: " << new_clus.energy() << endmsg;
      }
    }

    // That's all!

    return StatusCode::SUCCESS;
  }

  // get a map of mcID --> cluster
  // input: cluster_collections --> list of handles to all cluster collections
  std::map<eic::Index, eic::ConstCluster> indexedClusters(const eic::ClusterCollection& clusters) const {
    std::map<eic::Index, eic::ConstCluster> matched = {};
    for (const auto& cluster : clusters) {
      if (msgLevel(MSG::VERBOSE)) {
        verbose() << " --> Found cluster: " << cluster.id() << " with mcID " << cluster.mcID() << " and energy "
                  << cluster.energy() << endmsg;
      }
      if (!cluster.mcID()) {
        if (msgLevel(MSG::VERBOSE)) {
          verbose() << "   --> WARNING: no valid MC truth link found, skipping cluster..." << endmsg;
        }
        continue;
      }
      const bool duplicate = matched.count(cluster.mcID());
      if (duplicate) {
        if (msgLevel(MSG::VERBOSE)) {
          verbose() << "   --> WARNING: this is a duplicate mcID, keeping the higher energy cluster" << endmsg;
        }
        if (cluster.energy() < matched[cluster.mcID()].energy()) {
          continue;
        }
      }
      matched[cluster.mcID()] = cluster;
    }
    return matched;
  }
};
DECLARE_COMPONENT(TruthEnergyPositionClusterMerger)

} // namespace Jug::Fast
