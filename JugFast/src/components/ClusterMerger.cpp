// TODO needs full rework to run off mc-cluster relations instead
#if 0
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
  // Output
  DataHandle<eicd::ClusterCollection> m_outputClusters{"OutputClusters", Gaudi::DataHandle::Writer, this};
public:
  ClusterMerger(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputClusters", m_inputClusters, "Input cluster collection");
    declareProperty("outputClusters", m_outputClusters, "Cluster collection with good energy precision");
  }

  StatusCode initialize() override { return StatusCode::SUCCESS; }

  StatusCode execute() override {
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Merging cluster that belong to the same primary particle" << endmsg;
    }
    // input
    const auto& split = *(m_inputClusters.get());
    // output
    auto& merged    = *(m_outputClusters.createAndPut());

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
    auto clusterMap = indexedClusterLists(split);

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
          debug() << "   --> Only a single cluster, energy: " << clus.energy()
                  << " for this particle, copying" << endmsg;
        }
        merged.push_back(clus.clone());
      } else {
        auto new_clus = merged.create();
        // calculate aggregate info
        float energy      = 0;
        float energyError = 0;
        float time        = 0;
        int nhits = 0;
        eicd::VectorXYZ position;
        for (const auto& clus : clusters) {
          if (msgLevel(MSG::DEBUG)) {
            debug() << "   --> Adding cluster with energy: " << clus.energy() << endmsg;
          }
          energy += clus.energy();
          energyError += clus.energyError() * clus.energyError();
          time += clus.time() * clus.energy();
          nhits += clus.nhits();
          position = position.add(clus.position().scale(energy));
          new_clus.addclusters(clus);
          for (const auto& hit : clus.hits()) {
            new_clus.addhits(hit);
          }
        }
        new_clus.energy(energy);
        new_clus.energyError(sqrt(energyError));
        new_clus.time(time / energy);
        new_clus.nhits(nhits);
        new_clus.position(position.scale(1/energy));
        if (msgLevel(MSG::DEBUG)) {
          debug() << "   --> Merged cluster with energy: " << new_clus.energy() << endmsg;
        }
      }
    }

    // That's all!

    return StatusCode::SUCCESS;
  }

  // get a map of mcID --> std::vector<cluster> for clusters that belong together
  std::map<eicd::Index, std::vector<eicd::ConstCluster>> indexedClusterLists(const eicd::ClusterCollection& clusters) const {
    std::map<eicd::Index, std::vector<eicd::ConstCluster>> matched = {};
    for (const auto& cluster : clusters) {
      if (msgLevel(MSG::VERBOSE)) {
        verbose() << " --> Found cluster with mcID " << cluster.mcID() << " and energy "
                  << cluster.energy() << endmsg;
      }
      if (!cluster.mcID()) {
        if (msgLevel(MSG::VERBOSE)) {
          verbose() << "   --> WARNING: no valid MC truth link found, skipping cluster..." << endmsg;
        }
        continue;
      }
      if (!matched.count(cluster.mcID())) {
        matched[cluster.mcID()] = {};
      }
      matched[cluster.mcID()].push_back(cluster);
    }
    return matched;
  }
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ClusterMerger)

} // namespace Jug::Fast
#endif
