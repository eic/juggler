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
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/ReconstructedParticleRelationsCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/VectorPolar.h"

namespace Jug::Fast {

class MatchClusters : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  // input data
  DataHandle<dd4pod::Geant4ParticleCollection> m_inputMCParticles{"mcparticles", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_inputParticles{"ReconstructedChargedParticles",
                                                                    Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleRelationsCollection> m_inputRelations{"ReconstructedChargedParticlesRelations",
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
  DataHandle<eic::ReconstructedParticleRelationsCollection> m_outputRelations{"ReconstructedParticleRelations",
                                                                              Gaudi::DataHandle::Writer, this};

  MatchClusters(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputMCParticles", m_inputMCParticles, "mcparticles");
    declareProperty("inputParticles", m_inputParticles, "ReconstructedChargedParticles");
    declareProperty("inputRelations", m_inputRelations, "ReconstructedChargedParticleRelations");
    declareProperty("outputParticles", m_outputParticles, "ReconstructedParticles");
    declareProperty("outputRelations", m_outputRelations, "ReconstructedParticleRelations");
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
      debug() << "Adding cluster info to event" << endmsg;
    }
    // input collection
    const auto& mcparticles = *(m_inputMCParticles.get());
    const auto& inparts     = *(m_inputParticles.get());
    const auto& inrels      = *(m_inputRelations.get());
    auto& outparts          = *(m_outputParticles.createAndPut());
    auto& outrels           = *(m_outputRelations.createAndPut());

    // get an indexed map of all clusters
    auto ecalClusterMap = indexedClusters(m_inputEcalClusterCollections);
    auto hcalClusterMap = indexedClusters(m_inputHcalClusterCollections);

    // loop over all tracks and link matched clusters where applicable
    // (removing matched clusters from the cluster maps)
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Matching clusters to charged particles..." << endmsg;
    }
    for (size_t trkIdx = 0; trkIdx < inparts.size(); ++trkIdx) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "Processing charged particle " << trkIdx << endmsg;
      }
      outparts.push_back(inparts[trkIdx]);
      outrels.push_back(inrels[trkIdx]);
      auto rel = outrels[trkIdx];
      if (!rel.mcID()) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << " --> cannot match track without associated mcID" << endmsg;
        }
        continue;
      }
      if (ecalClusterMap.count(rel.mcID())) {
        const auto& clus = ecalClusterMap[rel.mcID()];
        if (msgLevel(MSG::DEBUG)) {
          debug() << " --> found matching Ecal cluster (" << clus.ID() << ")" << endmsg;
        }
        rel.ecalID(clus.ID());
        ecalClusterMap.erase(rel.mcID());
      }
      if (hcalClusterMap.count(rel.mcID())) {
        const auto& clus = hcalClusterMap[rel.mcID()];
        if (msgLevel(MSG::DEBUG)) {
          debug() << " --> found matching Hcal cluster (" << clus.ID() << ")" << endmsg;
        }
        rel.hcalID(clus.ID());
        hcalClusterMap.erase(rel.mcID());
      }
    }
    // Now loop over all remaining Ecal clusters and add neutrals. Also add in Hcal energy
    // if a matching cluster is available
    if (msgLevel(MSG::DEBUG)) {
      debug() << "Creating neutrals for remaining Ecal clusters..." << endmsg;
    }
    for (const auto& [mcID, eclus] : ecalClusterMap) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "Processing unmatched cluster (" << eclus.ID() << ")" << endmsg;
      }
      // get mass/PID from mcparticles, photon in case the matched particle is charged.
      const auto& mc    = mcparticles[mcID.value];
      const double mass = (!mc.charge()) ? mc.mass() : 0;
      const int32_t pid = (!mc.charge()) ? mc.pdgID() : 22;
      if (msgLevel(MSG::DEBUG)) {
        if (mc.charge()) {
          debug() << " --> associated mcparticle is not a neutral (PID: " << mc.pdgID()
                  << ", setting the reconstructed particle ID to photon" << endmsg;
        }
        debug() << " --> found matching associated mcparticle with PID: " << pid << endmsg;
      }
      // Do we also have any other matching clusters we need to add energy from (just
      // considering HCAL for now)
      const eic::ConstCluster* optional_hclus = nullptr;
      if (hcalClusterMap.count(mcID)) {
        optional_hclus = &(hcalClusterMap[mcID]);
      }
      // Reconstruct our neutrals and add them to the list
      const auto [part, rel] =
          reconstruct_neutral(eclus, optional_hclus, mass, pid, {static_cast<int32_t>(outparts.size()), algorithmID()});
      outparts.push_back(part);
      outrels.push_back(rel);
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
        matched[cluster.ID()] = cluster;
      }
    }
    return matched;
  }

  // reconstruct a neutral cluster
  // (for now assuming the vertex is at (0,0,0))
  std::pair<eic::ReconstructedParticle, eic::ReconstructedParticleRelations>
  reconstruct_neutral(const eic::ConstCluster& eclus, const eic::ConstCluster* optional_hclus, const double mass,
                      const int32_t pid, const eic::Index& newID) const {
    const float energy   = (optional_hclus) ? eclus.energy() + optional_hclus->energy() : eclus.energy();
    const float momentum = sqrt(energy * energy - mass * mass);
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
    // And relational info
    eic::ReconstructedParticleRelations rel;
    rel.recID(newID);
    rel.mcID(eclus.mcID());
    rel.ecalID(eclus.ID());
    if (optional_hclus) {
      rel.hcalID(optional_hclus->ID());
    }
    return {part, rel};
  }
};

DECLARE_COMPONENT(MatchClusters)

} // namespace Jug::Fast
