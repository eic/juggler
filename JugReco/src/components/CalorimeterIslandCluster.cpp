/*
 *  Island Clustering Algorithm for Calorimeter Blocks
 *  1. group all the adjacent modules
 *  2. split the groups between their local maxima with the energy deposit above <minClusterCenterEdep>
 *  Output hits collection with clusterID
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
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ProtoClusterCollection.h"
#include "eicd/VectorXY.h"
#include "eicd/VectorXYZ.h"

using namespace Gaudi::Units;

namespace Jug::Reco {
  // helper functions to get distance between hits
  static eic::VectorXY localDistXY(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return {h1.local().x - h2.local().x, h1.local().y - h2.local().y}; 
  }
  static eic::VectorXY localDistXZ(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return {h1.local().x - h2.local().x, h1.local().z - h2.local().z};
  }
  static eic::VectorXY localDistYZ(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return {h1.local().y - h2.local().y, h1.local().z - h2.local().z};
  }
  static eic::VectorXY dimScaledLocalDistXY(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return {2.*(h1.local().x - h2.local().x)/(h1.dimension().x + h2.dimension().x),
            2.*(h1.local().y - h2.local().y)/(h1.dimension().y + h2.dimension().y)};
  }
  static eic::VectorXY globalDistRPhi(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return {h1.position().r() - h2.position().r(), h1.position().phi() - h2.position().phi()};
  }
  static eic::VectorXY globalDistEtaPhi(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return {h1.position().eta() - h2.position().eta(), h1.position().phi() - h2.position().phi()};
  }
  // name: {method, units}
  static std::map<std::string, 
                  std::tuple<
                    std::function<eic::VectorXY(eic::ConstCalorimeterHit, eic::ConstCalorimeterHit)>,
                    std::vector<double>>> distMethods {
      {"localDistXY", {localDistXY, {mm, mm}}},
      {"localDistXZ", {localDistXZ, {mm, mm}}},
      {"localDistYZ", {localDistYZ, {mm, mm}}},
      {"dimScaledLocalDistXY", {dimScaledLocalDistXY, {1., 1.}}},
      {"globalDistRPhi", {globalDistRPhi, {mm, rad}}},
      {"globalDistEtaPhi", {globalDistEtaPhi, {1., rad}}},
  };

  /**
   *  Island Clustering Algorithm for Calorimeter Blocks.
   *
   *  1. group all the adjacent modules
   *  2. split the groups between their local maxima with the energy deposit above <minClusterCenterEdep>
   *  3. Output hits collection with clusterID
   *
   *  References:
   *      https://cds.cern.ch/record/687345/files/note01_034.pdf
   *      https://www.jlab.org/primex/weekly_meetings/primexII/slides_2012_01_20/island_algorithm.pdf
   *
   * \ingroup reco
   */
  class CalorimeterIslandCluster : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:
    Gaudi::Property<bool>                       m_splitCluster{this, "splitCluster", true};
    Gaudi::Property<double>                     m_minClusterHitEdep{this, "minClusterHitEdep", 0.};
    Gaudi::Property<double>                     m_minClusterCenterEdep{this, "minClusterCenterEdep", 50.0 * MeV};
    DataHandle<eic::CalorimeterHitCollection>   m_inputHitCollection{"inputHitCollection",
                                                                      Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ProtoClusterCollection>     m_outputProtoCollection{"outputProtoClusterCollection",
                                                                      Gaudi::DataHandle::Writer, this};

    // neighbour checking distances
    Gaudi::Property<double>                     m_sectorDist{this, "sectorDist", 5.0 * cm};
    Gaudi::Property<std::vector<double>>        u_localDistXY{this, "localDistXY", {}};
    Gaudi::Property<std::vector<double>>        u_localDistXZ{this, "localDistXZ", {}};
    Gaudi::Property<std::vector<double>>        u_localDistYZ{this, "localDistYZ", {}};
    Gaudi::Property<std::vector<double>>        u_globalDistRPhi{this, "globalDistRPhi", {}};
    Gaudi::Property<std::vector<double>>        u_globalDistEtaPhi{this, "globalDistEtaPhi", {}};
    Gaudi::Property<std::vector<double>>        u_dimScaledLocalDistXY{this, "dimScaledLocalDistXY", {1.8, 1.8}};
    // neighbor checking function
    std::function<eic::VectorXY(eic::ConstCalorimeterHit, eic::ConstCalorimeterHit)> hitsDist;

    // unitless counterparts of the input parameters
    double minClusterHitEdep, minClusterCenterEdep, sectorDist;
    std::array<double, 2> neighbourDist = {0., 0.};

    CalorimeterIslandCluster(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
        , AlgorithmIDMixin(name, info())
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputProtoClusterCollection", m_outputProtoCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }

      // unitless conversion, keep consistency with juggler internal units (GeV, mm, ns, rad)
      minClusterHitEdep = m_minClusterHitEdep.value() / GeV;
      minClusterCenterEdep = m_minClusterCenterEdep.value() / GeV;
      sectorDist = m_sectorDist.value() / mm;

      // set coordinate system
      auto set_dist_method = [this] (const Gaudi::Property<std::vector<double>> &uprop) {
        if (not uprop.size()) { return false; }
        auto &[method, units] = distMethods[uprop.name()];
        if (uprop.size() != units.size()) {
          info() << units.size() << endmsg;
          warning() << fmt::format("Expect {} values from {}, received {}: ({}), ignored it.", units.size(),
                                   uprop.name(), uprop.size(), fmt::join(uprop.value(), ", ")) << endmsg;
          return false;
        } else {
          for (size_t i = 0; i < units.size(); ++i) { neighbourDist[i] = uprop.value()[i]/units[i]; }
          hitsDist = method;
          info() << fmt::format("Clustering uses {} with distances <= [{}]",
                                uprop.name(), fmt::join(neighbourDist, ","))
                 << endmsg;
        }
        return true;
      };

      std::vector<Gaudi::Property<std::vector<double>>> uprops {
          u_localDistXY,
          u_localDistXZ,
          u_localDistYZ,
          u_globalDistRPhi,
          u_globalDistEtaPhi,
          // default one should be the last one
          u_dimScaledLocalDistXY,
      };

      bool method_found = false;
      for (auto &uprop : uprops) {
        if (set_dist_method(uprop)) { method_found = true; break; }
      }
      if (not method_found) {
        error() << "Cannot determine the clustering coordinates" << endmsg;
        return StatusCode::FAILURE;
      }

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& hits = *m_inputHitCollection.get();
      // Create output collections
      auto& proto = *m_outputProtoCollection.createAndPut();

      // group neighboring hits
      std::vector<std::vector<std::pair<uint32_t, eic::ConstCalorimeterHit>>> groups;

      std::vector<bool> visits(hits.size(), false);
      for (size_t i = 0; i < hits.size(); ++i) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << fmt::format("hit {:d}: energy = {:.4f} MeV, local = ({:.4f}, {:.4f}) mm, "
                                 "global=({:.4f}, {:.4f}, {:.4f}) mm, layer = {:d}, sector = {:d}.",
                                 i, hits[i].energy() * 1000., hits[i].local().x, hits[i].local().y,
                                 hits[i].position().x, hits[i].position().y, hits[i].position().z, hits[i].layer(),
                                 hits[i].sector())
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

      size_t clusterID = 0;
      for (auto& group : groups) {
        if (group.empty()) {
          continue;
        }
        auto maxima = find_maxima(group, !m_splitCluster.value());
        split_group(group, maxima, clusterID, proto);
        if (msgLevel(MSG::DEBUG)) {
          debug() << "hits in a group: " << group.size() << ", "
                  << "local maxima: " << maxima.size() << endmsg;
          debug() << "total number of clusters so far: " << clusterID << ", " << endmsg;
        }
      }

      return StatusCode::SUCCESS;
    }

  private:
    // helper function to group hits
    inline bool is_neighbour(const eic::ConstCalorimeterHit& h1, const eic::ConstCalorimeterHit& h2) const
    {
      // in the same sector
      if (h1.sector() == h2.sector()) {
        auto dist = hitsDist(h1, h2);
        return (dist.x <= neighbourDist[0]) && (dist.y <= neighbourDist[1]);
      // different sector, local coordinates do not work, using global coordinates
      } else {
        // sector may have rotation (barrel), so z is included
        return (h1.position().subtract(h2.position())).mag() <= sectorDist;
      }
    }

    // grouping function with Depth-First Search
    void dfs_group(std::vector <std::pair<uint32_t, eic::ConstCalorimeterHit>> & group, int idx,
                   const eic::CalorimeterHitCollection& hits, std::vector<bool>& visits) const {
      // not a qualified hit to particpate clustering, stop here
      if (hits[idx].energy() < minClusterHitEdep) {
        visits[idx] = true;
        return;
      }

      group.push_back({idx, hits[idx]});
      visits[idx] = true;
      for (size_t i = 0; i < hits.size(); ++i) {
        if (visits[i] || !is_neighbour(hits[idx], hits[i])) {
          continue;
        }
        dfs_group(group, i, hits, visits);
      }
    }

    // find local maxima that above a certain threshold
    std::vector<eic::ConstCalorimeterHit>
    find_maxima(const std::vector<std::pair<uint32_t, eic::ConstCalorimeterHit>>& group, bool global = false) const {
      std::vector<eic::ConstCalorimeterHit> maxima;
      if (group.empty()) {
        return maxima;
      }

      if (global) {
        int mpos = 0;
        for (size_t i = 0; i < group.size(); ++i) {
          if (group[mpos].second.energy() < group[i].second.energy()) {
            mpos = i;
          }
        }
        if (group[mpos].second.energy() >= minClusterCenterEdep) {
            maxima.push_back(group[mpos].second);
        }
        return maxima;
      }

      for (auto& [idx, hit] : group) {
        // not a qualified center
        if (hit.energy() < minClusterCenterEdep) {
          continue;
        }

        bool maximum = true;
        for (auto& [idx2, hit2] : group) {
          if (hit == hit2) {
            continue;
          }

          if (is_neighbour(hit, hit2) && hit2.energy() > hit.energy()) {
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
    inline void vec_normalize(std::vector<double>& vals) const
    {
      double total = 0.;
      for (auto& val : vals) {
        total += val;
      }
      for (auto& val : vals) {
        val /= total;
      }
    }

    // split a group of hits according to the local maxima
    void split_group(std::vector<std::pair<uint32_t, eic::ConstCalorimeterHit>>& group,
                     const std::vector<eic::ConstCalorimeterHit>& maxima, size_t& clusterID,
                     eic::ProtoClusterCollection& proto) const {
      // special cases
      if (maxima.size() == 0) {
        return;
      } else if (maxima.size() == 1) {
        eic::ProtoCluster pcl{{static_cast<int32_t>(clusterID), algorithmID()}};
        for (auto& [idx, hit] : group) {
          pcl.addhits({hit.ID(), idx, 1.});
        }
        clusterID += 1;
        return;
      }

      // split between maxima
      // TODO, here we can implement iterations with profile, or even ML for better splits
      std::vector<double> weights(maxima.size(), 1.);
      std::vector<eic::ProtoCluster> pcls;
      const size_t n_clus = clusterID + 1;
      for (size_t k = 0; k < maxima.size(); ++k) {
        pcls.push_back({{static_cast<int32_t>(n_clus + k), algorithmID()}});
      }

      size_t i      = 0;
      for (const auto& [idx, hit] : group) {
        size_t j     = 0;
        // calculate weights for local maxima
        for (const auto& chit : maxima) {
          double dist_ref = chit.dimension().x;
          double energy   = chit.energy();
          double dist     = hitsDist(chit, hit).mag();
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
          pcls[k].addhits({hit.ID(), idx, weight});
        }
        i += 1;
      }
      for (auto& pcl : pcls) {
        proto.push_back(pcl);
      }
      clusterID += maxima.size();
      return;
    }
  };

  DECLARE_COMPONENT(CalorimeterIslandCluster)

} // namespace Jug::Reco

