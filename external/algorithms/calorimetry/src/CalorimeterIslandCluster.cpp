// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Chao, Chao Peng, Wouter Deconinck, Jihee Kim, Whitney Armstrong

/*
 *  Island Clustering Algorithm for Calorimeter Blocks
 *  1. group all the adjacent modules
 *  2. split the groups between their local maxima with the energy deposit above <minClusterCenterEdep>
 *
 *  Author: Chao Peng (ANL), 09/27/2020
 *  References:
 *      https://cds.cern.ch/record/687345/files/note01_034.pdf
 *      https://www.jlab.org/primex/weekly_meetings/primexII/slides_2012_01_20/island_algorithm.pdf
 */
#include <algorithm>
#include <functional>
#include <tuple>

#include "fmt/format.h"

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

// Event Model related classes
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/ProtoClusterCollection.h"
#include "edm4eic/vector_utils.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Vector2f.h"

#if defined __has_include
#  if __has_include ("edm4eic/Vector3f.h")
#    include "edm4eic/Vector3f.h"
#  endif
#  if __has_include ("edm4eic/Vector2f.h")
#    include "edm4eic/Vector2f.h"
#  endif
#endif

namespace edm4eic {
  class Vector2f;
  class Vector3f;
}

using namespace Gaudi::Units;

namespace {

using CaloHit = edm4eic::CalorimeterHit;
using CaloHitCollection = edm4eic::CalorimeterHitCollection;

using Vector2f = std::conditional_t<
  std::is_same_v<decltype(edm4eic::CalorimeterHitData::position), edm4hep::Vector3f>,
  edm4hep::Vector2f,
  edm4eic::Vector2f
>;

// helper functions to get distance between hits
static Vector2f localDistXY(const CaloHit& h1, const CaloHit& h2) {
  const auto delta = h1.getLocal() - h2.getLocal();
  return {delta.x, delta.y};
}
static Vector2f localDistXZ(const CaloHit& h1, const CaloHit& h2) {
  const auto delta = h1.getLocal() - h2.getLocal();
  return {delta.x, delta.z};
}
static Vector2f localDistYZ(const CaloHit& h1, const CaloHit& h2) {
  const auto delta = h1.getLocal() - h2.getLocal();
  return {delta.y, delta.z};
}
static Vector2f dimScaledLocalDistXY(const CaloHit& h1, const CaloHit& h2) {
  const auto delta = h1.getLocal() - h2.getLocal();
  const auto dimsum = h1.getDimension() + h2.getDimension();
  return {2 * delta.x / dimsum.x, 2 * delta.y / dimsum.y};
}
static Vector2f globalDistRPhi(const CaloHit& h1, const CaloHit& h2) {
  using vector_type = decltype(Vector2f::a);
  return {
    static_cast<vector_type>(
      edm4eic::magnitude(h1.getPosition()) - edm4eic::magnitude(h2.getPosition())
    ),
    static_cast<vector_type>(
      edm4eic::angleAzimuthal(h1.getPosition()) - edm4eic::angleAzimuthal(h2.getPosition())
    )
  };
}
static Vector2f globalDistEtaPhi(const CaloHit& h1,
                                       const CaloHit& h2) {
  using vector_type = decltype(Vector2f::a);
  return {
    static_cast<vector_type>(
      edm4eic::eta(h1.getPosition()) - edm4eic::eta(h2.getPosition())
    ),
    static_cast<vector_type>(
      edm4eic::angleAzimuthal(h1.getPosition()) - edm4eic::angleAzimuthal(h2.getPosition())
    )
  };
}
// name: {method, units}
static std::map<std::string,
                std::tuple<std::function<Vector2f(const CaloHit&, const CaloHit&)>, std::vector<double>>>
    distMethods{
        {"localDistXY", {localDistXY, {mm, mm}}},        {"localDistXZ", {localDistXZ, {mm, mm}}},
        {"localDistYZ", {localDistYZ, {mm, mm}}},        {"dimScaledLocalDistXY", {dimScaledLocalDistXY, {1., 1.}}},
        {"globalDistRPhi", {globalDistRPhi, {mm, rad}}}, {"globalDistEtaPhi", {globalDistEtaPhi, {1., rad}}},
    };

} // namespace
namespace Jug::Reco {

/**
 *  Island Clustering Algorithm for Calorimeter Blocks.
 *
 *  1. group all the adjacent modules
 *  2. split the groups between their local maxima with the energy deposit above <minClusterCenterEdep>
 *
 *  References:
 *      https://cds.cern.ch/record/687345/files/note01_034.pdf
 *      https://www.jlab.org/primex/weekly_meetings/primexII/slides_2012_01_20/island_algorithm.pdf
 *
 * \ingroup reco
 */
class CalorimeterIslandCluster : public GaudiAlgorithm {
private:
  Gaudi::Property<bool> m_splitCluster{this, "splitCluster", true};
  Gaudi::Property<double> m_minClusterHitEdep{this, "minClusterHitEdep", 0.};
  Gaudi::Property<double> m_minClusterCenterEdep{this, "minClusterCenterEdep", 50.0 * MeV};
  DataHandle<CaloHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::ProtoClusterCollection> m_outputProtoCollection{"outputProtoClusterCollection",
                                                                  Gaudi::DataHandle::Writer, this};

  // neighbour checking distances
  Gaudi::Property<double> m_sectorDist{this, "sectorDist", 5.0 * cm};
  Gaudi::Property<std::vector<double>> u_localDistXY{this, "localDistXY", {}};
  Gaudi::Property<std::vector<double>> u_localDistXZ{this, "localDistXZ", {}};
  Gaudi::Property<std::vector<double>> u_localDistYZ{this, "localDistYZ", {}};
  Gaudi::Property<std::vector<double>> u_globalDistRPhi{this, "globalDistRPhi", {}};
  Gaudi::Property<std::vector<double>> u_globalDistEtaPhi{this, "globalDistEtaPhi", {}};
  Gaudi::Property<std::vector<double>> u_dimScaledLocalDistXY{this, "dimScaledLocalDistXY", {1.8, 1.8}};
  // neighbor checking function
  std::function<Vector2f(const CaloHit&, const CaloHit&)> hitsDist;

  // unitless counterparts of the input parameters
  double minClusterHitEdep{0}, minClusterCenterEdep{0}, sectorDist{0};
  std::array<double, 2> neighbourDist = {0., 0.};

public:
  CalorimeterIslandCluster(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputProtoClusterCollection", m_outputProtoCollection, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    // unitless conversion, keep consistency with juggler internal units (GeV, mm, ns, rad)
    minClusterHitEdep    = m_minClusterHitEdep.value() / GeV;
    minClusterCenterEdep = m_minClusterCenterEdep.value() / GeV;
    sectorDist           = m_sectorDist.value() / mm;

    // set coordinate system
    auto set_dist_method = [this](const Gaudi::Property<std::vector<double>>& uprop) {
      if (uprop.size() == 0) {
        return false;
      }
      auto& [method, units] = distMethods[uprop.name()];
      if (uprop.size() != units.size()) {
        info() << units.size() << endmsg;
        warning() << fmt::format("Expect {} values from {}, received {}: ({}), ignored it.", units.size(), uprop.name(),
                                 uprop.size(), fmt::join(uprop.value(), ", "))
                  << endmsg;
        return false;
      } else {
        for (size_t i = 0; i < units.size(); ++i) {
          neighbourDist[i] = uprop.value()[i] / units[i];
        }
        hitsDist = method;
        info() << fmt::format("Clustering uses {} with distances <= [{}]", uprop.name(), fmt::join(neighbourDist, ","))
               << endmsg;
      }
      return true;
    };

    std::vector<Gaudi::Property<std::vector<double>>> uprops{
        u_localDistXY,
        u_localDistXZ,
        u_localDistYZ,
        u_globalDistRPhi,
        u_globalDistEtaPhi,
        // default one should be the last one
        u_dimScaledLocalDistXY,
    };

    bool method_found = false;
    for (auto& uprop : uprops) {
      if (set_dist_method(uprop)) {
        method_found = true;
        break;
      }
    }
    if (not method_found) {
      error() << "Cannot determine the clustering coordinates" << endmsg;
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& hits = *(m_inputHitCollection.get());
    // Create output collections
    auto& proto = *(m_outputProtoCollection.createAndPut());

    // group neighboring hits
    std::vector<std::vector<std::pair<uint32_t, CaloHit>>> groups;

    std::vector<bool> visits(hits.size(), false);
    for (size_t i = 0; i < hits.size(); ++i) {
      if (msgLevel(MSG::DEBUG)) {
        const auto& hit = hits[i];
        debug() << fmt::format("hit {:d}: energy = {:.4f} MeV, local = ({:.4f}, {:.4f}) mm, "
                               "global=({:.4f}, {:.4f}, {:.4f}) mm",
                               i, hit.getEnergy() * 1000., hit.getLocal().x, hit.getLocal().y, hit.getPosition().x,
                               hit.getPosition().y, hit.getPosition().z)
                << endmsg;
      }
      // already in a group
      if (visits[i]) {
        continue;
      }
      groups.emplace_back();
      // create a new group, and group all the neighboring hits
      dfs_group(groups.back(), i, hits, visits);
    }

    for (auto& group : groups) {
      if (group.empty()) {
        continue;
      }
      auto maxima = find_maxima(group, !m_splitCluster.value());
      split_group(group, maxima, proto);
      if (msgLevel(MSG::DEBUG)) {
        debug() << "hits in a group: " << group.size() << ", "
                << "local maxima: " << maxima.size() << endmsg;
      }
    }

    return StatusCode::SUCCESS;
  }

private:
  // helper function to group hits
  inline bool is_neighbour(const CaloHit& h1, const CaloHit& h2) const {
    // in the same sector
    if (h1.getSector() == h2.getSector()) {
      auto dist = hitsDist(h1, h2);
      return (dist.a <= neighbourDist[0]) && (dist.b <= neighbourDist[1]);
      // different sector, local coordinates do not work, using global coordinates
    } else {
      // sector may have rotation (barrel), so z is included
      return (edm4eic::magnitude(h1.getPosition() - h2.getPosition()) <= sectorDist);
    }
  }

  // grouping function with Depth-First Search
  void dfs_group(std::vector<std::pair<uint32_t, CaloHit>>& group, int idx,
                 const CaloHitCollection& hits, std::vector<bool>& visits) const {
    // not a qualified hit to particpate clustering, stop here
    if (hits[idx].getEnergy() < minClusterHitEdep) {
      visits[idx] = true;
      return;
    }

    group.emplace_back(idx, hits[idx]);
    visits[idx] = true;
    for (size_t i = 0; i < hits.size(); ++i) {
      if (visits[i] || !is_neighbour(hits[idx], hits[i])) {
        continue;
      }
      dfs_group(group, i, hits, visits);
    }
  }

  // find local maxima that above a certain threshold
  std::vector<CaloHit>
  find_maxima(const std::vector<std::pair<uint32_t, CaloHit>>& group,
              bool global = false) const {
    std::vector<CaloHit> maxima;
    if (group.empty()) {
      return maxima;
    }

    if (global) {
      int mpos = 0;
      for (size_t i = 0; i < group.size(); ++i) {
        if (group[mpos].second.getEnergy() < group[i].second.getEnergy()) {
          mpos = i;
        }
      }
      if (group[mpos].second.getEnergy() >= minClusterCenterEdep) {
        maxima.push_back(group[mpos].second);
      }
      return maxima;
    }

    for (const auto& [idx, hit] : group) {
      // not a qualified center
      if (hit.getEnergy() < minClusterCenterEdep) {
        continue;
      }

      bool maximum = true;
      for (const auto& [idx2, hit2] : group) {
        if (hit == hit2) {
          continue;
        }

        if (is_neighbour(hit, hit2) && hit2.getEnergy() > hit.getEnergy()) {
          maximum = false;
          break;
        }
      }

      if (maximum) {
        maxima.push_back(hit);
      }
    }

    return maxima;
  }

  // helper function
  inline static void vec_normalize(std::vector<double>& vals) {
    double total = 0.;
    for (auto& val : vals) {
      total += val;
    }
    for (auto& val : vals) {
      val /= total;
    }
  }

  // split a group of hits according to the local maxima
  void split_group(std::vector<std::pair<uint32_t, CaloHit>>& group, const std::vector<CaloHit>& maxima,
                   edm4eic::ProtoClusterCollection& proto) const {
    // special cases
    if (maxima.empty()) {
      if (msgLevel(MSG::VERBOSE)) {
        verbose() << "No maxima found, not building any clusters" << endmsg;
      }
      return;
    } else if (maxima.size() == 1) {
      edm4eic::MutableProtoCluster pcl;
      for (auto& [idx, hit] : group) {
        pcl.addToHits(hit);
        pcl.addToWeights(1.);
      }
      proto.push_back(pcl);
      if (msgLevel(MSG::VERBOSE)) {
        verbose() << "A single maximum found, added one ProtoCluster" << endmsg;
      }
      return;
    }

    // split between maxima
    // TODO, here we can implement iterations with profile, or even ML for better splits
    std::vector<double> weights(maxima.size(), 1.);
    std::vector<edm4eic::MutableProtoCluster> pcls;
    for (size_t k = 0; k < maxima.size(); ++k) {
      pcls.emplace_back();
    }

    size_t i = 0;
    for (const auto& [idx, hit] : group) {
      size_t j = 0;
      // calculate weights for local maxima
      for (const auto& chit : maxima) {
        double dist_ref = chit.getDimension().x;
        double energy   = chit.getEnergy();
        double dist     = edm4eic::magnitude(hitsDist(chit, hit));
        weights[j]      = std::exp(-dist / dist_ref) * energy;
        j += 1;
      }

      // normalize weights
      vec_normalize(weights);

      // ignore small weights
      for (auto& w : weights) {
        if (w < 0.02) {
          w = 0;
        }
      }
      vec_normalize(weights);

      // split energy between local maxima
      for (size_t k = 0; k < maxima.size(); ++k) {
        double weight = weights[k];
        if (weight <= 1e-6) {
          continue;
        }
        pcls[k].addToHits(hit);
        pcls[k].addToWeights(weight);
      }
      i += 1;
    }
    for (auto& pcl : pcls) {
      proto.push_back(pcl);
    }
    if (msgLevel(MSG::VERBOSE)) {
      verbose() << "Multiple (" << maxima.size() << ") maxima found, added a ProtoClusters for each maximum" << endmsg;
    }
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(CalorimeterIslandCluster)

} // namespace Jug::Reco
