/*
 *  Reconstruct the cluster with Center of Gravity method
 *  Logarithmic weighting is used for mimicing energy deposit in transverse direction
 *
 *  Author: Chao Peng (ANL), 09/27/2020
 */
#include <algorithm>
#include <cstring>
#include <functional>
#include <map>

#include "boost/algorithm/string/join.hpp"
#include "fmt/format.h"
#include <boost/range/adaptor/map.hpp>

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
#include "eicd/ProtoClusterCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/Cluster2DInfoCollection.h"
#include "eicd/VectorPolar.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  // weighting functions (with place holders for hit energy, total energy, one parameter and module
  // type enum
  static double constWeight(double /*E*/, double /*tE*/, double /*p*/, int /*type*/) { return 1.0; }
  static double linearWeight(double E, double /*tE*/, double /*p*/, int /*type*/) { return E; }
  static double logWeight(double E, double tE, double base, int /*type*/)
  {
    return std::max(0., base + std::log(E / tE));
  }

  static const std::map<std::string, std::function<double(double, double, double, int)>>
    weightMethods {
      {"none", constWeight},
      {"linear", linearWeight},
      {"log", logWeight},
    };

  /** Clustering with center of gravity method.
   *
   *  Reconstruct the cluster with Center of Gravity method
   *  Logarithmic weighting is used for mimicing energy deposit in transverse direction
   *
   * \ingroup reco
   */
  class ClusterRecoCoG : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:
    Gaudi::Property<double>                   m_sampFrac{this, "samplingFraction", 1.0};
    Gaudi::Property<double>                   m_logWeightBase{this, "logWeightBase", 3.6};
    Gaudi::Property<double>                   m_depthCorrection{this, "depthCorrection", 0.0};
    Gaudi::Property<std::string>              m_energyWeight{this, "energyWeight", "log"};
    Gaudi::Property<std::string>              m_moduleDimZName{this, "moduleDimZName", ""};
    DataHandle<eic::CalorimeterHitCollection> m_inputHits{"inputHitCollection",
                                                          Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ProtoClusterCollection>   m_inputProto{"inputProtoClusterCollection",
                                                           Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ClusterCollection>        m_outputClusters{"outputClusterCollection",
                                                               Gaudi::DataHandle::Writer, this};
    DataHandle<eic::Cluster2DInfoCollection>  m_outputInfo{"outputInfoCollection",
                                                                  Gaudi::DataHandle::Writer, this};
    // Pointer to the geometry service
    SmartIF<IGeoSvc>                                   m_geoSvc;
    double                                             m_depthCorr;
    std::function<double(double, double, double, int)> weightFunc;

    ClusterRecoCoG(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin<>(name, info())
    {
      declareProperty("inputHitCollection", m_inputHits, "");
      declareProperty("inputProtoClusterCollection", m_inputProto, "");
      declareProperty("outputClusterCollection", m_outputClusters, "");
      declareProperty("outputInfoCollection", m_outputInfo, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }
      // update depth correction if a name is provided
      if (!m_moduleDimZName.value().empty()) {
        m_depthCorrection = m_geoSvc->detector()->constantAsDouble(m_moduleDimZName);
      }

      // select weighting method
      std::string ew = m_energyWeight.value();
      // make it case-insensitive
      std::transform(ew.begin(), ew.end(), ew.begin(), [](char s) { return std::tolower(s); });
      auto it = weightMethods.find(ew);
      if (it == weightMethods.end()) {
        error() << fmt::format(
                       "Cannot find energy weighting method {}, choose one from [{}]",
                       m_energyWeight,
                       boost::algorithm::join(weightMethods | boost::adaptors::map_keys, ", "))
                << endmsg;
        return StatusCode::FAILURE;
      }
      weightFunc = it->second;
      // info() << "z_length " << depth << endmsg;
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& hits     = *m_inputHits.get();
      const auto& proto    = *m_inputProto.get();
      auto&       clusters = *m_outputClusters.createAndPut();
      auto&       info     = *m_outputInfo.createAndPut();

      // Create a map of clusterID --> associated hits by looping over our clustered hits
      std::map<int, std::vector<std::pair<eic::ConstProtoCluster, 
                                          eic::ConstCalorimeterHit>>> cluster_map;
      for (const auto& pc : proto) {
        const size_t clusterID = pc.clusterID().value;
        if (!cluster_map.count(clusterID)) {
          cluster_map[clusterID] = {};
        }
        size_t idx;
        for (idx = 0; idx < hits.size(); ++idx) {
          if (hits[idx].ID() == pc.hitID()) {
            break;
          }
        }
        if (idx >= hits.size()) {
          continue;
        }
        cluster_map[clusterID].push_back({pc, hits[idx]});
      }

      for (const auto& [idx, hit_info] : cluster_map) {
        auto cl = reconstruct(hit_info, idx);

        if (msgLevel(MSG::DEBUG)) {
          debug() << cl.nhits() << " hits: " << cl.energy() / GeV << " GeV, (" << cl.position().x / mm << ", "
                  << cl.position().y / mm << ", " << cl.position().z / mm << ")" << endmsg;
        }
        clusters.push_back(cl);
        info.push_back({cl.ID(), cl.position(), cl.position().eta()});
      }

      return StatusCode::SUCCESS;
    }

  private:
    template <typename T1>
    eic::VectorPolar cart_to_polar(const T1& cart)
    {
      auto r = std::sqrt(cart.x * cart.x + cart.y * cart.y + cart.z * cart.z);
      return eic::VectorPolar{r, std::acos(cart.z / r), std::atan2(cart.y, cart.x)};
    }

    eic::Cluster reconstruct(const std::vector<std::pair<eic::ConstProtoCluster,
                                                         eic::ConstCalorimeterHit>>& hit_info, 
                             const int idx) const
    {
      eic::Cluster cl;
      cl.ID({idx, algorithmID()});
      cl.nhits(hit_info.size());

      // no hits
      if (msgLevel(MSG::DEBUG)) {
        debug() << "hit size = " << hit_info.size() << endmsg;
      }
      if (hit_info.empty()) {
        return cl;
      }

      // calculate total energy, find the cell with the maximum energy deposit
      float totalE   = 0.;
      float maxE     = 0.;
      auto time = hit_info[0].second.time();
      for (const auto& [proto, hit] : hit_info) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "hit energy = " << hit.energy() << " hit weight: " << proto.weight() << endmsg;
        }
        auto energy = hit.energy() * proto.weight();
        totalE += energy;
        if (energy > maxE) {
          maxE     = energy;
          time = hit.time();
        }
      }
      cl.energy(totalE/m_sampFrac);
      cl.energyError(0.);
      cl.time(time);

      // center of gravity with logarithmic weighting
      float tw = 0.;
      eic::VectorXYZ v;
      for (const auto& [proto, hit] : hit_info) {
        float w = weightFunc(hit.energy()*proto.weight(), totalE, m_logWeightBase.value(), 0);
        tw += w;
        v = v.add(hit.position().scale(w));
      }
      if (tw == 0.) {
        warning() << "zero total weights encountered, you may want to adjust your weighting parameter." << endmsg;
      }
      cl.position(v.scale(1/tw));
      cl.positionError({}); // @TODO: Covariance matrix

      return cl;
    }
  };

  DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

