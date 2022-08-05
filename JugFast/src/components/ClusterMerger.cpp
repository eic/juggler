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
#include "eicd/ClusterCollection.h"
#include "eicd/MCRecoClusterParticleAssociationCollection.h"
#include "eicd/vector_utils.h"

using namespace Gaudi::Units;

namespace Jug::Fast {

/** Simple algorithm to merge clusters orinating from the same particle together,
 *  based on the MC truth.
 *
 * \ingroup fast
 */
class ClusterMerger : public GaudiAlgorithm {
private:
  // Input
  DataHandle<eicd::ClusterCollection> m_inputClusters{"InputClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::MCRecoClusterParticleAssociationCollection> m_inputAssociations{"InputAssociations", Gaudi::DataHandle::Reader, this};
  // Output
  DataHandle<eicd::ClusterCollection> m_outputClusters{"OutputClusters", Gaudi::DataHandle::Writer, this};
  DataHandle<eicd::MCRecoClusterParticleAssociationCollection> m_outputAssociations{"OutputAssociations", Gaudi::DataHandle::Writer, this};
public:
  ClusterMerger(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputClusters", m_inputClusters, "Input cluster collection");
    declareProperty("inputAssociations", m_inputAssociations, "Input cluster association");
    declareProperty("outputClusters", m_outputClusters, "Cluster collection with good energy precision");
    declareProperty("outputAssociations", m_outputAssociations, "Cluster associations with good energy precision");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Merging cluster that belong to the same primary particle" << endmsg;
    }
    // input
    const auto& split = *(m_inputClusters.get());
    const auto& assoc = *(m_inputAssociations.get());
    // output
    auto& merged = *(m_outputClusters.createAndPut());
    auto& assoc2 = *(m_outputAssociations.createAndPut());

    if (!split.size()) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "Nothing to do for this event, returning..." << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 0/1: Getting indexed list of clusters..." << endmsg;
    }
    // get an indexed map of all vectors of clusters, indexed by mcID
    auto clusterMap = indexedClusterLists(split, assoc);

    // loop over all position clusters and match with energy clusters
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 1/1: Merging clusters where needed" << endmsg;
    }
    for (const auto& [mcID, clusters] : clusterMap) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Processing " << clusters.size() << " clusters for mcID " << mcID << endmsg;
      }
      if (clusters.size() == 1) {
        const auto& clus = clusters[0];
        if (msgLevel(MSG::DEBUG)) {
          debug() << "   --> Only a single cluster, energy: " << clus.getEnergy()
                  << " for this particle, copying" << endmsg;
        }
        auto new_clus = clus.clone();
        merged.push_back(new_clus);
        auto ca = assoc2.create();
        ca.setRecID(new_clus.getObjectID().index);
        ca.setSimID(mcID);
        ca.setWeight(1.0);
        ca.setRec(new_clus);
        //ca.setSim(//FIXME);
      } else {
        auto new_clus = merged.create();
        // calculate aggregate info
        float energy      = 0;
        float energyError = 0;
        float time        = 0;
        int nhits = 0;
        eicd::Vector3f position;
        for (const auto& clus : clusters) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "   --> Adding cluster with energy: " << clus.getEnergy() << endmsg;
          }
          energy += clus.getEnergy();
          energyError += clus.getEnergyError() * clus.getEnergyError();
          time += clus.getTime() * clus.getEnergy();
          nhits += clus.getNhits();
          position = position + energy * clus.getPosition();
          new_clus.addToClusters(clus);
          for (const auto& hit : clus.getHits()) {
            new_clus.addToHits(hit);
          }
        }
        new_clus.setEnergy(energy);
        new_clus.setEnergyError(sqrt(energyError));
        new_clus.setTime(time / energy);
        new_clus.setNhits(nhits);
        new_clus.setPosition(position / energy);
        if (msgLevel(MSG::DEBUG)) {
          debug() << "   --> Merged cluster with energy: " << new_clus.getEnergy() << endmsg;
        }
        auto ca = assoc2.create();
        ca.setSimID(mcID);
        ca.setWeight(1.0);
        ca.setRec(new_clus);
      }
    }

    // That's all!

    return StatusCode::SUCCESS;
  }

  // get a map of MCParticle index--> std::vector<Cluster> for clusters that belong together
  std::map<int, std::vector<eicd::Cluster>> indexedClusterLists(
      const eicd::ClusterCollection& clusters,
      const eicd::MCRecoClusterParticleAssociationCollection& associations
  ) const {

    std::map<int, std::vector<eicd::Cluster>> matched = {};

    // loop over clusters
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
        verbose() << " --> Found cluster with mcID " << mcID << " and energy "
                  << cluster.getEnergy() << endmsg;
      }

      if (mcID < 0) {
        if (msgLevel(MSG::VERBOSE)) {
          verbose() << "   --> WARNING: no valid MC truth link found, skipping cluster..." << endmsg;
        }
        continue;
      }

      if (!matched.count(mcID)) {
        matched[mcID] = {};
      }
      matched[mcID].push_back(cluster);
    }
    return matched;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ClusterMerger)

} // namespace Jug::Fast
