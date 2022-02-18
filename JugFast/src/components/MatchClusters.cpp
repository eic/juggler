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
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/VectorPolar.h"

namespace Jug::Fast {

class MatchClusters : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  // input data
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"MCParticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_inputParticles{"ReconstructedChargedParticles",
                                                                    Gaudi::DataHandle::Reader, this};
  Gaudi::Property<std::vector<std::string>> m_inputEcalClusters{
      this, "inputEcalClusters", {}, "Ecal clusters to be aggregated"};
  Gaudi::Property<std::vector<std::string>> m_inputHcalClusters{
      this, "inputHcalClusters", {}, "Hcal clusters to be aggregated"};
  std::vector<DataHandle<eic::ClusterCollection>*> m_inputEcalClusterCollections;
  std::vector<DataHandle<eic::ClusterCollection>*> m_inputHcalClusterCollections;

  // output data
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"ReconstructedParticles",
                                                                     Gaudi::DataHandle::Writer, this};

  MatchClusters(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputMCParticles, "MCParticles");
    declareProperty("inputParticles", m_inputParticles, "ReconstructedChargedParticles");
    declareProperty("outputParticles", m_outputParticles, "ReconstructedParticles");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    m_inputEcalClusterCollections = getClusterCollections(m_inputEcalClusters);
    m_inputHcalClusterCollections = getClusterCollections(m_inputHcalClusters);
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
    auto ecalClusterMap = indexedClusters(m_inputEcalClusterCollections);
    auto hcalClusterMap = indexedClusters(m_inputHcalClusterCollections);

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
      if (ecalClusterMap.count(part.mcID())) {
        const auto& clus = ecalClusterMap[part.mcID()];
        if (msgLevel(MSG::DEBUG)) {
          debug() << "    --> found matching Ecal cluster (" << clus.ID() << "), energy: " << clus.energy() << endmsg;
        }
        part.ecalID(clus.ID());
        ecalClusterMap.erase(part.mcID());
      }
      if (hcalClusterMap.count(part.mcID())) {
        const auto& clus = hcalClusterMap[part.mcID()];
        if (msgLevel(MSG::DEBUG)) {
          debug() << " --> found matching Hcal cluster (" << clus.ID() << "), energy: " << clus.energy() << endmsg;
        }
        part.hcalID(clus.ID());
        hcalClusterMap.erase(part.mcID());
      }
    }
    // Now loop over all remaining Ecal clusters and add neutrals. Also add in Hcal energy
    // if a matching cluster is available
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Step 2/2: Creating neutrals for remaining Ecal clusters..." << endmsg;
    }
    for (const auto& [mcID, eclus] : ecalClusterMap) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << " --> Processing unmatched ECAL cluster (" << eclus.ID() << "), energy: " << eclus.energy()
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
      // Do we also have any other matching clusters we need to add energy from (just
      // considering HCAL for now)
      const eic::ConstCluster* optional_hclus = nullptr;
      if (hcalClusterMap.count(mcID)) {
        optional_hclus = &(hcalClusterMap[mcID]);
        if (msgLevel(MSG::DEBUG)) {
          debug() << "   --> found matching HCAL cluster (" << *optional_hclus
                  << "), energy: " << optional_hclus->energy() << endmsg;
        }
      }
      // Reconstruct our neutrals and add them to the list
      const auto part =
          reconstruct_neutral(eclus, optional_hclus, mass, pid, {static_cast<int32_t>(outparts.size()), algorithmID()});
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
          verbose() << " --> Found cluster: " << cluster.ID() << " with mcID " << cluster.mcID() << " and energy "
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
  eic::ReconstructedParticle reconstruct_neutral(const eic::ConstCluster& eclus,
                                                 const eic::ConstCluster* optional_hclus, const double mass,
                                                 const int32_t pid, const eic::Index& newID) const {
    const float energy   = (optional_hclus) ? eclus.energy() + optional_hclus->energy() : eclus.energy();
    const float momentum = energy < mass ? 0 : std::hypot(energy, mass);
    const eic::VectorPolar p{momentum, eclus.position().theta(), eclus.position().phi()};
    // setup our particle
    eic::ReconstructedParticle part;
    part.ID(newID);
    part.p(p);
    part.time(eclus.time());
    part.pid(pid);
    part.status(0);
    part.charge(0);
    part.direction({p.theta, p.phi});
    part.momentum(momentum);
    part.energy(energy);
    part.mass(mass);
    part.mcID(eclus.mcID());
    part.ecalID(eclus.ID());
    if (optional_hclus) {
      part.hcalID(optional_hclus->ID());
    }
    return part;
  }
}; // namespace Jug::Fast

DECLARE_COMPONENT(MatchClusters)

} // namespace Jug::Fast
