// TODO needs full rework to run off list of mc-cluster relations instead
#if 0
// Takes a list of particles (presumed to be from tracking), and all available clusters.
// 1. Match clusters to their tracks using the mcID field
// 2. For unmatched clusters create neutrals and add to the particle list

#include <algorithm>
#include <cmath>

#include <fmt/format.h>

#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/vector_utils.h"

namespace Jug::Fast {

class MatchClusters : public GaudiAlgorithm {
public:
  // input data
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"MCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_inputParticles{"ReconstructedChargedParticles",
                                                                    Gaudi::DataHandle::Reader, this};
  Gaudi::Property<std::vector<std::string>> m_inputClusters{this, "inputClusters", {}, "Clusters to be aggregated"};
  std::vector<DataHandle<eic::ClusterCollection>*> m_inputClusterCollections;

  // output data
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"ReconstructedParticles",
                                                                     Gaudi::DataHandle::Writer, this};

  MatchClusters(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticles, "MCParticles");
    declareProperty("inputParticles", m_inputParticles, "ReconstructedChargedParticles");
    declareProperty("outputParticles", m_outputParticles, "ReconstructedParticles");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    m_inputClusterCollections = getClusterCollections(m_inputClusters);
    return StatusCode::SUCCESS;
  }
  StatusCode execute() override {
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Processing cluster info for new event" << endmsg;
    }
    // input collection
    const auto& mcparticles = *(m_inputMCParticles.get());
    const auto& inparts     = *(m_inputParticles.get());
    auto& outparts          = *(m_outputParticles.createAndPut());

    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 0/2: Getting indexed list of clusters..." << endmsg;
    }
    // get an indexed map of all clusters
    auto clusterMap = indexedClusters(m_inputClusterCollections);

    // loop over all tracks and link matched clusters where applicable
    // (removing matched clusters from the cluster maps)
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 1/2: Matching clusters to charged particles..." << endmsg;
    }
    for (size_t trkIdx = 0; trkIdx < inparts.size(); ++trkIdx) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Processing charged particle " << trkIdx << ", PID: " << inparts[trkIdx].pid()
                << ", energy: " << inparts[trkIdx].energy() << endmsg;
      }
      outparts.push_back(inparts[trkIdx].clone());
      auto part = outparts[trkIdx];
      if (!part.mcID()) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "    --> cannot match track without associated mcID" << endmsg;
        }
        continue;
      }
      if (clusterMap.count(part.mcID())) {
        const auto& clus = clusterMap[part.mcID()];
        if (msgLevel(MSG::DEBUG)) {
          debug() << "    --> found matching cluster with energy: " << clus.energy() << endmsg;
        }
        clusterMap.erase(part.mcID());
      }
    }
    // Now loop over all remaining clusters and add neutrals. Also add in Hcal energy
    // if a matching cluster is available
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 2/2: Creating neutrals for remaining clusters..." << endmsg;
    }
    for (const auto& [mcID, clus] : clusterMap) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Processing unmatched cluster with energy: " << clus.energy()
                << endmsg;
      }

      // get mass/PID from mcparticles, 0 (unidentified) in case the matched particle is charged.
      const auto& mc    = mcparticles[mcID.value];
      const double mass = (!mc.getCharge()) ? mc.getMass() : 0;
      const int32_t pid = (!mc.getCharge()) ? mc.getPDG() : 0;
      if (msgLevel(MSG::DEBUG)) {
        if (mc.getCharge()) {
          debug() << "   --> associated mcparticle is not a neutral (PID: " << mc.getPDG()
                  << "), setting the reconstructed particle ID to 0 (unidentified)" << endmsg;
        }
        debug() << "   --> found matching associated mcparticle with PID: " << pid << ", energy: " << mc.getEnergy()
                << endmsg;
      }
      // Reconstruct our neutrals and add them to the list
      const auto part = reconstruct_neutral(clus, mass, pid);
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Reconstructed neutral particle with PID: " << part.pid() << ", energy: " << part.energy()
                << endmsg;
      }
      outparts.push_back(part);
      // That's all!
    }
    return StatusCode::SUCCESS;
  }

private:
  std::vector<DataHandle<eic::ClusterCollection>*> getClusterCollections(const std::vector<std::string>& cols) {
    std::vector<DataHandle<eic::ClusterCollection>*> ret;
    for (const auto& colname : cols) {
      debug() << "initializing cluster collection: " << colname << endmsg;
      ret.push_back(new DataHandle<eic::ClusterCollection>{colname, Gaudi::DataHandle::Reader, this});
    }
    return ret;
  }

  // get a map of mcID --> cluster
  // input: cluster_collections --> list of handles to all cluster collections
  std::map<eic::Index, eic::ConstCluster>
  indexedClusters(const std::vector<DataHandle<eic::ClusterCollection>*>& cluster_collections) const {
    std::map<eic::Index, eic::ConstCluster> matched = {};
    for (const auto& cluster_handle : cluster_collections) {
      const auto& clusters = *(cluster_handle->get());
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
    }
    return matched;
  }

  // reconstruct a neutral cluster
  // (for now assuming the vertex is at (0,0,0))
  eic::ReconstructedParticle reconstruct_neutral(const eic::ConstCluster& clus, const double mass,
                                                 const int32_t pid) const {
    const float energy   = clus.energy();
    const float momentum = energy < mass ? 0 : std::sqrt(energy * energy - mass * mass);
    const auto p = eicd::normalizeVector(clus.position(), momentum);
    // setup our particle
    eic::ReconstructedParticle part;
    part.p(p);
    part.time(clus.time());
    part.pid(pid);
    part.status(0);
    part.charge(0);
    part.theta(p.theta());
    part.phi(p.phi());
    part.momentum(momentum);
    part.energy(energy);
    part.mass(mass);
    part.mcID(clus.mcID());
    return part;
  }
}; // namespace Jug::Fast

DECLARE_COMPONENT(MatchClusters)

} // namespace Jug::Fast
#endif
