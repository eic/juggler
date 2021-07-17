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

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"

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

  class ClusterRecoCoG : public GaudiAlgorithm {
  public:
    Gaudi::Property<double>                   m_sampFrac{this, "samplingFraction", 1.0};
    Gaudi::Property<double>                   m_logWeightBase{this, "logWeightBase", 3.6};
    Gaudi::Property<double>                   m_depthCorrection{this, "depthCorrection", 0.0};
    Gaudi::Property<std::string>              m_energyWeight{this, "energyWeight", "log"};
    Gaudi::Property<std::string>              m_moduleDimZName{this, "moduleDimZName", ""};
    DataHandle<eic::CalorimeterHitCollection> m_inputHits{"inputHitCollection",
                                                          Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ClusterCollection>        m_outputClusters{"outputClusterCollection",
                                                        Gaudi::DataHandle::Writer, this};
    // Pointer to the geometry service
    SmartIF<IGeoSvc>                                   m_geoSvc;
    double                                             m_depthCorr;
    std::function<double(double, double, double, int)> weightFunc;

    ClusterRecoCoG(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHits, "");
      declareProperty("outputClusterCollection", m_outputClusters, "");
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
      auto&       clusters = *m_outputClusters.createAndPut();

      // Create a map of clusterID --> associated hits by looping over our clustered hits
      std::map<int, std::vector<eic::ConstCalorimeterHit>> cluster_map;
      for (const auto& hit : hits) {
        const size_t idx = hit.clusterID();
        if (!cluster_map.count(idx)) {
          cluster_map[idx] = {};
        }
        cluster_map[idx].push_back(hit);
      }

      for (const auto& [idx, clhits] : cluster_map) {
        auto         clhit = reconstruct(clhits);
        eic::Cluster cl{idx,
                        static_cast<float>(clhit.energy() / m_sampFrac),
                        clhit.energy(),
                        static_cast<int>(clhits.size()),
                        clhit.position(),
                        {},
                        cart_to_polar(clhit.position()),
                        0.,
                        0.};

        debug() << cl.nhits() << " hits: " << cl.energy() / GeV << " GeV, (" << cl.position().x / mm
                << ", " << cl.position().y / mm << ", " << cl.position().z / mm << ")" << endmsg;
        clusters.push_back(cl);
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

    eic::CalorimeterHit reconstruct(const std::vector<eic::ConstCalorimeterHit>& hits) const
    {
      eic::CalorimeterHit res;
      // no hits
      debug() << "hit size = " << hits.size() << endmsg;
      if (hits.empty()) {
        return res;
        ;
      }

      // calculate total energy, find the cell with the maximum energy deposit
      float totalE   = 0.;
      float maxE     = 0.;
      auto  centerID = hits.begin()->cellID();
      for (const auto& hit : hits) {
        debug() << "hit energy = " << hit.energy() << endmsg;
        auto energy = hit.energy();
        totalE += energy;
        if (energy > maxE) {
          maxE     = energy;
          centerID = hit.cellID();
        }
      }
      res.cellID(centerID);
      res.energy(totalE);

      // center of gravity with logarithmic weighting
      float tw = 0., x = 0., y = 0., z = 0.;
      for (auto& hit : hits) {
        // suppress low energy contributions
        // info() << std::log(hit.energy()/totalE) << endmsg;
        float w = weightFunc(hit.energy(), totalE, m_logWeightBase.value(), 0);
        tw += w;
        x += hit.x() * w;
        y += hit.y() * w;
        z += hit.z() * w;
        /*
        debug() << hit.cellID() << ": (" << hit.local_x() << ", " << hit.local_y() << ", "
                << hit.local_z() << "), "
                << "(" << hit.x() << ", " << hit.y() << ", " << hit.z() << "), " << endmsg;
        */
      }
      if (tw == 0.) {
        warning() << "zero total weights encountered, you may want to adjust your weighting parameter." << endmsg;
      }
      res.position({x / tw, y / tw, z / tw});
      // convert global position to local position, use the cell with max edep as a reference
      const auto volman    = m_geoSvc->detector()->volumeManager();
      const auto alignment = volman.lookupDetElement(centerID).nominal();
      const auto lpos      = alignment.worldToLocal(dd4hep::Position(res.x(), res.y(), res.z()));

      // TODO: may need convert back to have depthCorrection in global positions
      res.local({lpos.x(), lpos.y(), lpos.z() + m_depthCorrection});
      return res;
    }
  };

  DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

