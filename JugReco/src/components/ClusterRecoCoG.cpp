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
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ProtoClusterCollection.h"
#include "eicd/ClusterCollection.h"
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
    // Constrain the cluster position eta to be within
    // the eta of the contributing hits. This is useful to avoid edge effects
    // for endcaps.
    Gaudi::Property<bool>                     m_enableEtaBounds{this, "enableEtaBounds", false};

    Gaudi::Property<std::string>              m_mcHits{this, "mcHits", ""};

    DataHandle<eic::CalorimeterHitCollection> m_inputHits{"inputHitCollection",
                                                          Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ProtoClusterCollection>   m_inputProto{"inputProtoClusterCollection",
                                                           Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ClusterCollection>        m_outputClusters{"outputClusterCollection",
                                                               Gaudi::DataHandle::Writer, this};

    // Monte Carlo particle source identifier
    const int32_t m_kMonteCarloSource{uniqueID<int32_t>("MCParticles")};
    // Optional handle to MC hits
    std::unique_ptr<DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputMC;

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
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      // Initialize the MC input hit collection if requested
      if (m_mcHits != "") {
        m_inputMC =
            std::make_unique<DataHandle<edm4hep::SimCalorimeterHitCollection>>(m_mcHits, Gaudi::DataHandle::Reader, this);
      }
      //
      //
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
      const auto& hits   = *m_inputHits.get();
      const auto& proto  = *m_inputProto.get();
      auto& clusters     = *m_outputClusters.createAndPut();
      // Optional MC data
      const edm4hep::SimCalorimeterHitCollection* mcHits = nullptr;
      if (m_inputMC) {
        mcHits = m_inputMC->get();
      }

      for (const auto& pcl : proto) {
        auto cl = reconstruct(pcl, hits, mcHits);

        if (msgLevel(MSG::DEBUG)) {
          debug() << cl.nhits() << " hits: " << cl.energy() / GeV << " GeV, (" << cl.position().x / mm << ", "
                  << cl.position().y / mm << ", " << cl.position().z / mm << ")" << endmsg;
        }
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

    eic::Cluster reconstruct(const eic::ConstProtoCluster& pcl, const eic::CalorimeterHitCollection& hits,
                             const edm4hep::SimCalorimeterHitCollection* /* mcHits */) const {
      eic::Cluster cl;
      cl.ID({pcl.ID(), algorithmID()});
      cl.nhits(pcl.hits_size());

      // no hits
      if (msgLevel(MSG::DEBUG)) {
        debug() << "hit size = " << pcl.hits_size() << endmsg;
      }
      if (pcl.hits_size() == 0) {
        return cl;
      }

      // calculate total energy, find the cell with the maximum energy deposit
      float totalE   = 0.;
      float maxE     = 0.;
      // Used to optionally constrain the cluster eta to those of the contributing hits
      float minHitEta = std::numeric_limits<float>::max();
      float maxHitEta = std::numeric_limits<float>::min();
      auto time = hits[pcl.hits(0).index].time();
      for (const auto& clhit : pcl.hits()) {
        const auto& hit = hits[clhit.index];
        if (msgLevel(MSG::DEBUG)) {
          debug() << "hit energy = " << hit.energy() << " hit weight: " << clhit.weight << endmsg;
        }
        auto energy = hit.energy() * clhit.weight;
        totalE += energy;
        if (energy > maxE) {
          maxE     = energy;
          time = hit.time();
        }
        const float eta = hit.position().eta();
        if (eta < minHitEta) {
          minHitEta = eta;
        }
        if (eta > maxHitEta) {
          maxHitEta = eta;
        }
      }
      cl.energy(totalE/m_sampFrac);
      cl.energyError(0.);
      cl.time(time);

      // center of gravity with logarithmic weighting
      float tw = 0.;
      eic::VectorXYZ v;
      for (const auto& clhit : pcl.hits()) {
        const auto& hit = hits[clhit.index];
        float w = weightFunc(hit.energy()*clhit.weight, totalE, m_logWeightBase.value(), 0);
        tw += w;
        v = v.add(hit.position().scale(w));
      }
      if (tw == 0.) {
        warning() << "zero total weights encountered, you may want to adjust your weighting parameter." << endmsg;
      }
      cl.position(v.scale(1/tw));
      cl.positionError({}); // @TODO: Covariance matrix

      // Optionally constrain the cluster to the hit eta values
      if (m_enableEtaBounds) {
        const bool overflow  = (cl.position().eta() > maxHitEta);
        const bool underflow = (cl.position().eta() < minHitEta);
        if (overflow || underflow) {
          const double newEta = overflow ? maxHitEta : minHitEta;
          const double newTheta = 2 * atan(exp(-newEta));
          const eic::VectorPolar oldPos = cl.position();
          cl.position(eic::VectorPolar(oldPos.r, newTheta, oldPos.phi));
          if (msgLevel(MSG::DEBUG)) {
            debug() << "Bound cluster position to contributing hits due to " << (overflow ? "overflow" : "underflow")
                    << endmsg;
          }
        }
      }

      // Additional convenience variables
      cl.polar(cl.position());
      cl.eta(cl.position().eta());

      // best estimate on the cluster direction is the cluster position
      // for simple 2D CoG clustering
      cl.direction({cl.polar().theta, cl.polar().phi});

      // Calculate radius
      // @TODO: add skewness
      if (cl.nhits() > 1) {
        double radius = 0;
        for (const auto& clhit : pcl.hits()) {
          const auto& hit  = hits[clhit.index];
          const auto delta = cl.position().subtract(hit.position());
          radius += delta.dot(delta);
        }
        radius = sqrt((1. / (cl.nhits() - 1.)) * radius);
        cl.radius(radius);
        // skewness not yet calculated
        cl.skewness(0);
      }

      // Optionally store the MC truth associated with the first hit in this cluster
      // FIXME no connection between cluster and truth in edm4hep
      //if (mcHits) {
      //  const auto& mc_hit    = (*mcHits)[pcl.hits(0).ID.value];
      //  cl.mcID({mc_hit.truth().trackID, m_kMonteCarloSource});
      //}

      return cl;
    }
  };

  DECLARE_COMPONENT(ClusterRecoCoG)

  } // namespace Jug::Reco

