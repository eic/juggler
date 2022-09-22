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
#include "fmt/ranges.h"

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

#ifdef USE_SYCL
    #include <CL/sycl.hpp>
    #include <unordered_map>
#endif

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

#ifdef USE_SYCL

// Functor throws async device errors to user thread
struct async_err_handler {
  void operator()(sycl::exception_list elist){
    for(auto e : elist){
      std::rethrow_exception(e);
    }
  }
};

// SYCL Initializations
sycl::default_selector device_selector;
sycl::queue queue(device_selector, async_err_handler{});

// enum to choose hitsDist method within SYCL code
enum class sycl_hits_dist {
    localXY,
    localXZ,
    localYZ,
    dimScaleLocalXY,
    globalRPhi,
    globalEtaPhi
};

sycl_hits_dist sHitsDist;
bool is_debug = false;

#endif

} // namespace
namespace Jug::Reco {

#ifdef USE_SYCL
  SYCL_EXTERNAL inline bool is_neighbour( const sycl::accessor<sycl::float3,1,sycl::access::mode::read>& lpos,
                                          const sycl::accessor<sycl::float3,1,sycl::access::mode::read>& gpos,
                                          const sycl::accessor<sycl::float3,1,sycl::access::mode::read>& dims,
                                          const sycl::accessor<int,1,sycl::access::mode::read>& sectors,
                                          const sycl::accessor<double,1,sycl::access::mode::read>& secDist,
                                          const sycl::accessor<double,1,sycl::access::mode::read>& neighDist,
                                          sycl::id<1> hit1, int hit2, int hitDist_type) {
    
    if(sectors[hit1] != sectors[hit2]) {
      auto delta = gpos[hit1] - gpos[hit2];
      return (sycl::length(delta) <= secDist[0]);
    } else {
      auto delta = lpos[hit1] - lpos[hit2];
      switch(hitDist_type){
        case 0: {
          // LocalXY
          return (delta.x() <= neighDist[0] && delta.y() <= neighDist[1]);
          break;
        }
        case 1: {
          // LocalXZ
          return (delta.x() <= neighDist[0] && delta.z() <= neighDist[1]);
          break;
        }
        case 2: {
          // LocalYZ
          return (delta.y() <= neighDist[0] && delta.z() <= neighDist[1]);
          break;
        }
        case 3: {
          // dimScaleLocalXY
          auto dimsum = dims[hit1] + dims[hit2];
          return ((2*delta.x()/dimsum.x()) <= neighDist[0] && (2*delta.y()/dimsum.y()) <= neighDist[1]);
          break;
        }
        case 4: {
          // globalRPhi
          double r = sycl::length(gpos[hit1]) - sycl::length(gpos[hit2]);
          double phi = sycl::atan2(gpos[hit1].y(), gpos[hit1].x()) - sycl::atan2(gpos[hit2].y(), gpos[hit2].x());
          return (r <= neighDist[0] && phi <= neighDist[1]);
          break;
        }
        case 5: {
          // globalEtaPhi
          double ang_polar = sycl::atan2(sycl::hypot(gpos[hit1].x(), gpos[hit1].y()), gpos[hit1].z());
          double eta = -sycl::log(sycl::tan(0.5 * ang_polar));
          double phi = sycl::atan2(gpos[hit1].y(), gpos[hit1].x()) - sycl::atan2(gpos[hit2].y(), gpos[hit2].x());
          return (eta <= neighDist[0] && phi <= neighDist[1]); 
        }
      }
    }

  }

  void parallel_group(std::vector<std::vector<std::pair<uint32_t, CaloHit>>>& groups,
                      const CaloHitCollection& hits, std::array<double,2>& neighbourDist,
                      double minClusterHitEdep, double sectorDist) {
    // Corner Cases
    if(hits.size() <= 0) return;

    // Host memory
    //
    // Get location data from hits
    std::vector<int32_t> sectors;
    std::vector<float> energy;
    std::vector<sycl::vec<float, 3>> lpos, gpos, dims;

    // Neighbour Indices
    std::vector<int> nidx (hits.size());

    // Can't filter out hits here as index numbers will not match
    for(size_t i = 0; i < hits.size(); i++) {
      lpos.emplace_back(hits[i].getLocal().x, hits[i].getLocal().y, hits[i].getLocal().z);
      gpos.emplace_back(hits[i].getPosition().x, hits[i].getPosition().y, hits[i].getPosition().z);
      dims.emplace_back(hits[i].getDimension().x, hits[i].getDimension().y, hits[i].getDimension().z);
      energy.emplace_back(hits[i].getEnergy());
      sectors.emplace_back(hits[i].getSector());
    }

    {
      // Device memory
      sycl::buffer<int, 1> nidx_buf (nidx.data(), sycl::range<1>(nidx.size()));
      sycl::buffer<double, 1> minClusterHitEdep_buf (&minClusterHitEdep, sycl::range<1>(1));
      sycl::buffer<float, 1> energy_buf (energy.data(), sycl::range<1>(energy.size()));

      sycl::buffer<sycl::float3,1> lpos_buf (lpos.data(), sycl::range<1>(lpos.size()));
      sycl::buffer<sycl::float3,1> gpos_buf (gpos.data(), sycl::range<1>(gpos.size()));
      sycl::buffer<sycl::float3,1> dims_buf (dims.data(), sycl::range<1>(dims.size()));

      sycl::buffer<int,1> sectors_buf (sectors.data(), sycl::range<1>(sectors.size()));
      sycl::buffer<double, 1> sectorDist_buf (&sectorDist, sycl::range<1>(1));
      sycl::buffer<double, 1> neighbourDist_buf (neighbourDist.data(), sycl::range<1>(neighbourDist.size()));

      try {
            // Initalize Neighbour Indices
            queue.submit([&](sycl::handler& h) {
              auto nidx_acc = nidx_buf.get_access<sycl::access::mode::write>(h);
              h.parallel_for(sycl::range<1>(hits.size()), [=](sycl::id<1> idx) {
                  nidx_acc[idx] = idx;
              });
            });

            /**
             * @brief Assign current vertex id (idx) to neighbour if idx < id
             * held at neighbour index. Emulates sequential assignment of 
             * clusters by DFS in parallel.
             */
            queue.submit([&](sycl::handler& h) {
              auto nidx_acc = nidx_buf.get_access<sycl::access::mode::atomic>(h);
              auto minClusterHitEdep_acc = minClusterHitEdep_buf.get_access<sycl::access::mode::read>(h);
              auto energy_acc = energy_buf.get_access<sycl::access::mode::read>(h);

              int hitsDist_type = static_cast<int>(sHitsDist);
              auto lpos_acc = lpos_buf.get_access<sycl::access::mode::read>(h);
              auto gpos_acc = gpos_buf.get_access<sycl::access::mode::read>(h);
              auto dims_acc = dims_buf.get_access<sycl::access::mode::read>(h);

              auto sectors_acc = sectors_buf.get_access<sycl::access::mode::read>(h);
              auto sectorDist_acc = sectorDist_buf.get_access<sycl::access::mode::read>(h);
              auto neighbourDist_acc = neighbourDist_buf.get_access<sycl::access::mode::read>(h);

              h.parallel_for(sycl::range<1>(hits.size()), [=](sycl::id<1> idx) {
                // not a qualified hit - skip
                if(energy_acc[idx] < minClusterHitEdep_acc[0]){
                  return;
                }

                for(size_t i = 0; i < nidx_acc.size(); i++) {

                  if(energy_acc[i] < minClusterHitEdep_acc[0]) continue;

                  if(!Jug::Reco::is_neighbour(lpos_acc,gpos_acc,dims_acc,sectors_acc,
                                   sectorDist_acc,neighbourDist_acc,idx,i,hitsDist_type)){
                    continue;
                  }
                  
                  // Atomic exchange of min element between current neighbour index nidx_acc[i] and idx
                  nidx_acc[i].fetch_min(idx);

                }

              });
            }).wait_and_throw();
      } catch(sycl::exception e) {
        std::cerr << "Caught SYCL Exception: " << e.what() << "\n";
      }

    } // Sync Device and Host memory

    if(is_debug){
        fmt::print("Parallel grouping results are:\n{}\n",fmt::join(nidx," "));
    }
    
    // Emplace index array into groups for further processing
    std::unordered_map<int, std::vector<std::pair<uint32_t, CaloHit>>> grp_map;
    for(size_t i = 0; i < hits.size(); i++) {
      grp_map[nidx[i]].emplace_back(i, hits[i]);
    }
    for(auto i : grp_map) {
      groups.emplace_back(i.second);
    }

  }

  std::vector<CaloHit>
  parallel_find_maxima(const std::vector<std::pair<uint32_t, CaloHit>>& group,
              std::array<double,2>& neighbourDist,
              double minClusterCenterEdep, double sectorDist, bool global = false) {
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

    // Prep Host memory for device offload
    std::vector<float> energy;
    std::vector<uint8_t> max_idx (group.size());

    // Get location data from hits
    std::vector<int32_t> sectors;
    std::vector<sycl::vec<float, 3>> lpos, gpos, dims;

    for (const auto& [idx, hit] : group){
      lpos.emplace_back(hit.getLocal().x, hit.getLocal().y, hit.getLocal().z);
      gpos.emplace_back(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
      dims.emplace_back(hit.getDimension().x, hit.getDimension().y, hit.getDimension().z);
      energy.push_back(hit.getEnergy());
      sectors.push_back(hit.getSector());
    }

    // Device memory
    {

      sycl::buffer<float, 1> energy_buf (energy.data(), sycl::range<1>(energy.size()));
      sycl::buffer<uint8_t, 1> max_idx_buf (max_idx.data(), sycl::range<1>(max_idx.size()));
      sycl::buffer<double, 1> minClusterCenterEdep_buf (&minClusterCenterEdep, sycl::range<1>(1));

      int hitsDist_type = static_cast<int>(sHitsDist);
      sycl::buffer<sycl::float3,1> lpos_buf (lpos.data(), sycl::range<1>(lpos.size()));
      sycl::buffer<sycl::float3,1> gpos_buf (gpos.data(), sycl::range<1>(gpos.size()));
      sycl::buffer<sycl::float3,1> dims_buf (dims.data(), sycl::range<1>(dims.size()));

      sycl::buffer<int32_t,1> sectors_buf (sectors.data(), sycl::range<1>(sectors.size()));
      sycl::buffer<double, 1> sectorDist_buf (&sectorDist, sycl::range<1>(1));
      sycl::buffer<double, 1> neighbourDist_buf (neighbourDist.data(), sycl::range<1>(neighbourDist.size()));
      
      try {
        queue.submit([&](sycl::handler& h){
          auto energy_acc = energy_buf.get_access<sycl::access::mode::read>(h);
          auto max_idx_acc = max_idx_buf.get_access<sycl::access::mode::write>(h);
          auto minClusterCenterEdep_acc = minClusterCenterEdep_buf.get_access<sycl::access::mode::read>(h);

          auto lpos_acc = lpos_buf.get_access<sycl::access::mode::read>(h);
          auto gpos_acc = gpos_buf.get_access<sycl::access::mode::read>(h);
          auto dims_acc = dims_buf.get_access<sycl::access::mode::read>(h);

          auto sectors_acc = sectors_buf.get_access<sycl::access::mode::read>(h);
          auto sectorDist_acc = sectorDist_buf.get_access<sycl::access::mode::read>(h);
          auto neighbourDist_acc = neighbourDist_buf.get_access<sycl::access::mode::read>(h);

          h.parallel_for(sycl::range<1>(group.size()), [=](sycl::id<1> idx){
            // not a qualified hit - skip
            if(energy_acc[idx] < minClusterCenterEdep_acc[0]){
              return;
            }
            
            bool is_max = true;

            for(size_t i = 0; i < energy_acc.size(); i++){
              if(idx == i){
                continue;
              }

              if(energy_acc[i] > energy_acc[idx]){
                if(Jug::Reco::is_neighbour(lpos_acc,gpos_acc,dims_acc,sectors_acc,
                                sectorDist_acc,neighbourDist_acc,idx,i,hitsDist_type)){
                      is_max = false;
                      break;
                  }
                      
              }
            }

            if(is_max){
              max_idx_acc[idx] = 1;
            }
            
          });
        }).wait_and_throw();

      } catch(sycl::exception e) {
        std::cerr << "Caught SYLC Exception: " << e.what() << "\n";
      }

    } // Sync Device memory to Host

    // Convert maxima index array to hit vector for further processing
    for(size_t i = 0; i < max_idx.size(); i++){
      if(max_idx[i] == 1){
        if(is_debug){
            std::cout << "Found Maxima at: " << group[i].first << "\n";
        }
        maxima.push_back(group[i].second);
      }
    }

    return maxima;
  }

#endif

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

        #ifdef USE_SYCL
        if(uprop.name() == "localDistXY")
            sHitsDist = sycl_hits_dist::localXY;
        if(uprop.name() == "localDistXZ")
            sHitsDist = sycl_hits_dist::localXZ;
        if(uprop.name() == "localDistYZ")
            sHitsDist = sycl_hits_dist::localYZ;
        if(uprop.name() == "dimScaledLocalDistXY")
            sHitsDist = sycl_hits_dist::dimScaleLocalXY;
        if(uprop.name() == "globalDistRPhi")
            sHitsDist = sycl_hits_dist::globalRPhi;
        if(uprop.name() == "globalDistEtaPhi")
            sHitsDist = sycl_hits_dist::globalEtaPhi;
        
        if(msgLevel(MSG::DEBUG)) is_debug = true;
        #endif

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

    #ifdef USE_SYCL
    for (size_t i = 0; i < hits.size(); ++i) {
      if (msgLevel(MSG::DEBUG)) {
        const auto& hit = hits[i];
        debug() << fmt::format("hit {:d}: energy = {:.4f} MeV, local = ({:.4f}, {:.4f}) mm, "
                               "global=({:.4f}, {:.4f}, {:.4f}) mm",
                               i, hit.getEnergy() * 1000., hit.getLocal().x, hit.getLocal().y, hit.getPosition().x,
                               hit.getPosition().y, hit.getPosition().z)
                << endmsg;
      }
    }

    parallel_group(groups, hits, neighbourDist, minClusterHitEdep, sectorDist);
    #else
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
    #endif

    for (auto& group : groups) {
      if (group.empty()) {
        continue;
      }
    #ifdef USE_SYCL
      auto maxima = parallel_find_maxima(group, neighbourDist, minClusterCenterEdep, sectorDist, !m_splitCluster.value());
    #else
      auto maxima = find_maxima(group, !m_splitCluster.value());
    #endif
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

  // find local maxima that are above a certain threshold
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
