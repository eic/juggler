// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten

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
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/MCRecoClusterParticleAssociationCollection.h"
#include <edm4eic/vector_utils.h>

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
private:
  // Input
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"MCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::ClusterCollection> m_energyClusters{"EnergyClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::MCRecoClusterParticleAssociationCollection> m_energyAssociations{"EnergyAssociations", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::ClusterCollection> m_positionClusters{"PositionClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::MCRecoClusterParticleAssociationCollection> m_positionAssociations{"PositionAssociations", Gaudi::DataHandle::Reader, this};
  // Output
  DataHandle<edm4eic::ClusterCollection> m_outputClusters{"OutputClusters", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4eic::MCRecoClusterParticleAssociationCollection> m_outputAssociations{"OutputAssociations", Gaudi::DataHandle::Writer, this};

public:
  TruthEnergyPositionClusterMerger(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticles, "MCParticles");
    declareProperty("inputEnergyClusters", m_energyClusters, "Cluster collection with good energy precision");
    declareProperty("inputEnergyAssociations", m_energyAssociations, "Cluster association with good energy precision");
    declareProperty("inputPositionClusters", m_positionClusters, "Cluster collection with good position precision");
    declareProperty("inputPositionAssociations", m_positionAssociations, "Cluster association with good position precision");
    declareProperty("outputClusters", m_outputClusters, "Cluster collection with good energy precision");
    declareProperty("outputAssociations", m_outputAssociations, "Cluster association with good energy precision");
  }

  StatusCode initialize() override { return StatusCode::SUCCESS; }

  StatusCode execute() override {
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Merging energy and position clusters for new event" << endmsg;
    }
    // input
    const auto& mcparticles  = *(m_inputMCParticles.get());
    const auto& energy_clus  = *(m_energyClusters.get());
    const auto& energy_assoc = *(m_energyAssociations.get());
    const auto& pos_clus     = *(m_positionClusters.get());
    const auto& pos_assoc    = *(m_positionAssociations.get());
    // output
    auto& merged_clus  = *(m_outputClusters.createAndPut());
    auto& merged_assoc = *(m_outputAssociations.createAndPut());

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
    auto energyMap = indexedClusters(energy_clus, energy_assoc);
    auto posMap    = indexedClusters(pos_clus, pos_assoc);

    // loop over all position clusters and match with energy clusters
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 1/2: Matching all position clusters to the available energy clusters..." << endmsg;
    }
    for (const auto& [mcID, pclus] : posMap) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Processing position cluster " << pclus.id() << ", mcID: " << mcID << ", energy: " << pclus.getEnergy()
                << endmsg;
      }
      if (energyMap.count(mcID)) {
        const auto& eclus = energyMap[mcID];
        auto new_clus = merged_clus.create();
        new_clus.setEnergy(eclus.getEnergy());
        new_clus.setEnergyError(eclus.getEnergyError());
        new_clus.setTime(pclus.getTime());
        new_clus.setNhits(pclus.getNhits() + eclus.getNhits());
        new_clus.setPosition(pclus.getPosition());
        new_clus.setPositionError(pclus.getPositionError());
        new_clus.addToClusters(pclus);
        new_clus.addToClusters(eclus);
        for (const auto& cl : {pclus, eclus}) {
          for (const auto& hit : cl.getHits()) {
            new_clus.addToHits(hit);
          }
          new_clus.addToSubdetectorEnergies(cl.getEnergy());
        }
        for (const auto& param : pclus.getShapeParameters()) {
          new_clus.addToShapeParameters(param);
        }
        if (msgLevel(MSG::DEBUG)) {
          debug() << "   --> Found matching energy cluster " << eclus.id() << ", energy: " << eclus.getEnergy() << endmsg;
          debug() << "   --> Created new combined cluster " << new_clus.id() << ", energy: " << new_clus.getEnergy() << endmsg;
        }

        // set association
        edm4eic::MutableMCRecoClusterParticleAssociation clusterassoc;
        clusterassoc.setRecID(new_clus.getObjectID().index);
        clusterassoc.setSimID(mcID);
        clusterassoc.setWeight(1.0);
        clusterassoc.setRec(new_clus);
        //clusterassoc.setSim(mcparticles[mcID]);
        merged_assoc.push_back(clusterassoc);

        // erase the energy cluster from the map, so we can in the end account for all
        // remaining clusters
        energyMap.erase(mcID);
      } else {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "   --> No matching energy cluster found, copying over position cluster" << endmsg;
        }
        auto new_clus = pclus.clone();
        new_clus.addToClusters(pclus);
        merged_clus.push_back(new_clus);

        // set association
        edm4eic::MutableMCRecoClusterParticleAssociation clusterassoc;
        clusterassoc.setRecID(new_clus.getObjectID().index);
        clusterassoc.setSimID(mcID);
        clusterassoc.setWeight(1.0);
        clusterassoc.setRec(new_clus);
        //clusterassoc.setSim(mcparticles[mcID]);
        merged_assoc.push_back(clusterassoc);
      }
    }
    // Collect remaining energy clusters. Use mc truth position for these clusters, as
    // they should really have a match in the position clusters (and if they don't it due
    // to a clustering error).
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 2/2: Collecting remaining energy clusters..." << endmsg;
    }
    for (const auto& [mcID, eclus] : energyMap) {
      const auto& mc = mcparticles[mcID];
      const auto& p = mc.getMomentum();
      const auto phi = std::atan2(p.y, p.x);
      const auto theta = std::atan2(std::hypot(p.x, p.y), p.z);
      auto new_clus = merged_clus.create();
      new_clus.setEnergy(eclus.getEnergy());
      new_clus.setEnergyError(eclus.getEnergyError());
      new_clus.setTime(eclus.getTime());
      new_clus.setNhits(eclus.getNhits());
      // use nominal radius of 110cm, and use start vertex theta and phi
      new_clus.setPosition(edm4eic::sphericalToVector(110.*cm, theta, phi));
      new_clus.addToClusters(eclus);
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Processing energy cluster " << eclus.id() << ", mcID: " << mcID << ", energy: " << eclus.getEnergy()
                << endmsg;
        debug() << "   --> Created new 'combined' cluster " << new_clus.id() << ", energy: " << new_clus.getEnergy() << endmsg;
      }

      // set association
      edm4eic::MutableMCRecoClusterParticleAssociation clusterassoc;
      clusterassoc.setRecID(new_clus.getObjectID().index);
      clusterassoc.setSimID(mcID);
      clusterassoc.setWeight(1.0);
      clusterassoc.setRec(new_clus);
      //clusterassoc.setSim(mc);
      merged_assoc.push_back(clusterassoc);
    }

    // That's all!
    return StatusCode::SUCCESS;
  }

  // get a map of MCParticle index --> cluster
  // input: cluster_collections --> list of handles to all cluster collections
  std::map<int, edm4eic::Cluster> indexedClusters(
      const edm4eic::ClusterCollection& clusters,
      const edm4eic::MCRecoClusterParticleAssociationCollection& associations
  ) const {

    std::map<int, edm4eic::Cluster> matched = {};

    for (const auto& cluster : clusters) {
      int mcID = -1;

      // find associated particle
      for (const auto& assoc : associations) {
        if (assoc.getRec() == cluster) {
          mcID = assoc.getSimID();
          break;
        }
      }

      if (msgLevel(MSG::VERBOSE)) {
        verbose() << " --> Found cluster: " << cluster.getObjectID().index << " with mcID " << mcID << " and energy "
                  << cluster.getEnergy() << endmsg;
      }

      if (mcID < 0) {
        if (msgLevel(MSG::VERBOSE)) {
          verbose() << "   --> WARNING: no valid MC truth link found, skipping cluster..." << endmsg;
        }
        continue;
      }

      const bool duplicate = matched.count(mcID);
      if (duplicate) {
        if (msgLevel(MSG::VERBOSE)) {
          verbose() << "   --> WARNING: this is a duplicate mcID, keeping the higher energy cluster" << endmsg;
        }
        if (cluster.getEnergy() < matched[mcID].getEnergy()) {
          continue;
        }
      }

      matched[mcID] = cluster;
    }
    return matched;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(TruthEnergyPositionClusterMerger)

} // namespace Jug::Fast
