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
#include "Math/Point2D.h"
#include "Math/Point3D.h"

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

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;
typedef ROOT::Math::XYPoint Point;
typedef ROOT::Math::XYZPoint Point3D;

namespace Jug::Reco {
  // helper functions to get distance between hits
  static Point localDistXY(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return Point(h1.local().local_x - h2.local().local_x, h1.local().local_y - h2.local().local_y);
  }
  static Point localDistXZ(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return Point(h1.local().local_x - h2.local().local_x, h1.local().local_z - h2.local().local_z);
  }
  static Point localDistYZ(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return Point(h1.local().local_y - h2.local().local_y, h1.local().local_z - h2.local().local_z);
  }
  static Point dimScaledLocalDistXY(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    return Point(2.*(h1.local().local_x - h2.local().local_x)/(h1.dimension().dim_x + h2.dimension().dim_x),
                 2.*(h1.local().local_y - h2.local().local_y)/(h1.dimension().dim_y + h2.dimension().dim_y));
  }
  static Point globalDistRPhi(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    Point3D p1(h1.position().x, h1.position().y, h1.position().z), p2(h2.position().x, h2.position().y, h2.position().z);
    return Point(p1.r() - p2.r(), p1.phi() - p2.phi());
  }
  static Point globalDistEtaPhi(eic::ConstCalorimeterHit h1, eic::ConstCalorimeterHit h2) {
    Point3D p1(h1.position().x, h1.position().y, h1.position().z), p2(h2.position().x, h2.position().y, h2.position().z);
    return Point(p1.eta() - p2.eta(), p1.phi() - p2.phi());
  }
  // name: {method, units}
  static std::map< std::string, std::tuple<
    std::function<Point(eic::ConstCalorimeterHit, eic::ConstCalorimeterHit)>,
    std::vector<double> > > distMethods {
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
  class CalorimeterIslandCluster : public GaudiAlgorithm {
  public:
    Gaudi::Property<bool>                       m_splitCluster{this, "splitCluster", true};
    Gaudi::Property<double>                     m_minClusterHitEdep{this, "minClusterHitEdep", 0.};
    Gaudi::Property<double>                     m_minClusterCenterEdep{this, "minClusterCenterEdep", 50.0 * MeV};
    DataHandle<eic::CalorimeterHitCollection>   m_inputHitCollection{"inputHitCollection",
                                                                      Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection>   m_splitHitCollection{"outputHitCollection",
                                                                      Gaudi::DataHandle::Writer, this};

    // neighbour checking distances
    Gaudi::Property<double>                     m_sectorDist{this, "sectorDist", 5.0 * cm};
    Gaudi::Property<std::vector<double>>        u_localDistXY{this, "localDistXY", {}};
    Gaudi::Property<std::vector<double>>        u_localDistXZ{this, "localDistXZ", {}};
    Gaudi::Property<std::vector<double>>        u_localDistYZ{this, "localDistYZ", {}};
    Gaudi::Property<std::vector<double>>        u_globalDistRPhi{this, "globalDistRPhi", {}};
    Gaudi::Property<std::vector<double>>        u_globalDistEtaPhi{this, "globalDistEtaPhi", {}};
    Gaudi::Property<std::vector<double>>        u_dimScaledLocalDistXY{this, "dimScaledLocalDistXY", {1.8, 1.8}};
    // neightbor checking function
    std::function<Point(eic::ConstCalorimeterHit, eic::ConstCalorimeterHit)> hitsDist;

    // unitless counterparts of the input parameters
    double minClusterHitEdep, minClusterCenterEdep, sectorDist;
    std::array<double, 2> neighbourDist = {0., 0.};

    CalorimeterIslandCluster(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_splitHitCollection, "");
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
      auto& clustered_hits = *m_splitHitCollection.createAndPut();

      // group neighboring hits
      std::vector<std::vector<eic::CalorimeterHit>> groups;

      std::vector<bool> visits(hits.size(), false);
      for (size_t i = 0; i < hits.size(); ++i) {
        debug() << fmt::format("hit {:d}: energy = {:.4f} MeV, local = ({:.4f}, {:.4f}) mm, "
                               "global=({:.4f}, {:.4f}, {:.4f}) mm, layer = {:d}, sector = {:d}.",
                               i, hits[i].energy()*1000., hits[i].local().local_x, hits[i].local().local_y,
                               hits[i].position().x, hits[i].position().y, hits[i].position().z,
                               hits[i].layerID(), hits[i].sectorID()) << endmsg;
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
        split_group(group, maxima, clusterID, clustered_hits);
        debug() << "hits in a group: " << group.size() << ", "
                << "local maxima: " << maxima.size() << endmsg;
        debug() << "total number of clusters so far: " << clusterID << ", " << endmsg;
      }

      return StatusCode::SUCCESS;
    }

  private:
    template <typename T>
    inline T dist3d(T t1, T t2, T t3) const
    {
      return std::sqrt(t1 * t1 + t2 * t2 + t3 * t3);
    }
    // helper function to group hits
    inline bool is_neighbour(const eic::ConstCalorimeterHit& h1, const eic::ConstCalorimeterHit& h2) const
    {
      // in the same sector
      if (h1.sectorID() == h2.sectorID()) {
        auto dist = hitsDist(h1, h2);
        return (dist.x() <= neighbourDist[0]) && (dist.y() <= neighbourDist[1]);
      // different sector, local coordinates do not work, using global coordinates
      } else {
        // sector may have rotation (barrel), so z is included
        return dist3d(h1.position().x - h2.position().x, h1.position().y - h2.position().y, h1.position().z - h2.position().z) <= sectorDist;
      }
    }

    // grouping function with Depth-First Search
    void dfs_group(std::vector<eic::CalorimeterHit>& group, int idx,
                   const eic::CalorimeterHitCollection& hits, std::vector<bool>& visits) const
    {
      // not a qualified hit to particpate clustering, stop here
      if (hits[idx].energy() < minClusterHitEdep) {
        visits[idx] = true;
        return;
      }

      eic::CalorimeterHit hit{hits[idx].cellID(),    hits[idx].clusterID(),
                              hits[idx].layerID(),   hits[idx].sectorID(),
                              hits[idx].energy(),    hits[idx].time(),
                              hits[idx].position(),  hits[idx].local(),
                              hits[idx].dimension(), 1};
      group.push_back(hit);
      visits[idx] = true;
      for (size_t i = 0; i < hits.size(); ++i) {
        if (visits[i] || !is_neighbour(hit, hits[i])) {
          continue;
        }
        dfs_group(group, i, hits, visits);
      }
    }

    // find local maxima that above a certain threshold
    std::vector<eic::ConstCalorimeterHit> find_maxima(const std::vector<eic::CalorimeterHit>& group,
                                                      bool global = false) const
    {
      std::vector<eic::ConstCalorimeterHit> maxima;
      if (group.empty()) {
        return maxima;
      }

      if (global) {
        int mpos = 0;
        for (size_t i = 0; i < group.size(); ++i) {
          if (group[mpos].energy() < group[i].energy()) {
            mpos = i;
          }
        }
        if (group[mpos].energy() >= minClusterCenterEdep) {
            maxima.push_back(group[mpos]);
        }
        return maxima;
      }

      for (auto& hit : group) {
        // not a qualified center
        if (hit.energy() < minClusterCenterEdep) {
          continue;
        }

        bool maximum = true;
        for (auto& hit2 : group) {
          if (hit == hit2)
            continue;

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
    // split_hits is used to persistify the data
    void split_group(std::vector<eic::CalorimeterHit>& group,
                     const std::vector<eic::ConstCalorimeterHit>& maxima,
                     size_t& clusterID, eic::CalorimeterHitCollection& clustered_hits) const
    {
      // special cases
      if (maxima.size() == 0) {
        return;
      } else if (maxima.size() == 1) {
        for (auto& hit : group) {
          hit.clusterID(clusterID);
          clustered_hits.push_back(hit);
        }
        clusterID += 1;
        return;
      }

      // split between maxima
      // TODO, here we can implement iterations with profile, or even ML for better splits
      std::vector<double>       weights(maxima.size());
      std::vector<eic::Cluster> splits(maxima.size());
      size_t                    n_clus = clusterID + 1;
      size_t                    i      = 0;
      for (auto it = group.begin(); it != group.end(); ++it, ++i) {
        auto   hedep = it->energy();
        size_t j     = 0;
        // calculate weights for local maxima
        for (auto cit = maxima.begin(); cit != maxima.end(); ++cit, ++j) {
          double dist_ref = cit->dimension().dim_x;
          double energy   = cit->energy();
          double dist     = hitsDist(*cit, *it).r();
          weights[j]      = std::exp(-dist / dist_ref) * energy;
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
        for (size_t k = 0; k < weights.size(); ++k) {
          double weight = weights[k];
          if (weight <= 1e-6) {
            continue;
          }
          auto hit = it->clone();
          hit.energy(hedep * weight);
          hit.type(1);
          hit.clusterID(n_clus + k);
          clustered_hits.push_back(hit);
        }
      }
      clusterID += splits.size();
      return;
    }
  };

  DECLARE_COMPONENT(CalorimeterIslandCluster)

} // namespace Jug::Reco

