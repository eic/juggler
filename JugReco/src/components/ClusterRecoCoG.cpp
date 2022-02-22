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

// Event Model related classes
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ProtoClusterCollection.h"
#include "eicd/vector_utils.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

// weighting functions (with place holders for hit energy, total energy, one parameter and module
// type enum
static double constWeight(double /*E*/, double /*tE*/, double /*p*/, int /*type*/) { return 1.0; }
static double linearWeight(double E, double /*tE*/, double /*p*/, int /*type*/) { return E; }
static double logWeight(double E, double tE, double base, int /*type*/) {
  return std::max(0., base + std::log(E / tE));
}

static const std::map<std::string, std::function<double(double, double, double, int)>> weightMethods{
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
class ClusterRecoCoG : public GaudiAlgorithm {
public:
  Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};
  Gaudi::Property<double> m_logWeightBase{this, "logWeightBase", 3.6};
  Gaudi::Property<double> m_depthCorrection{this, "depthCorrection", 0.0};
  Gaudi::Property<std::string> m_energyWeight{this, "energyWeight", "log"};
  Gaudi::Property<std::string> m_moduleDimZName{this, "moduleDimZName", ""};
  // Constrain the cluster position eta to be within
  // the eta of the contributing hits. This is useful to avoid edge effects
  // for endcaps.
  Gaudi::Property<bool> m_enableEtaBounds{this, "enableEtaBounds", false};

  //DataHandle<eic::CalorimeterHitCollection> m_inputHits{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ProtoClusterCollection> m_inputProto{"inputProtoClusterCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ClusterCollection> m_outputClusters{"outputClusterCollection", Gaudi::DataHandle::Writer, this};

  // Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  double m_depthCorr;
  std::function<double(double, double, double, int)> weightFunc;

  ClusterRecoCoG(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    //declareProperty("inputHitCollection", m_inputHits, "");
    declareProperty("inputProtoClusterCollection", m_inputProto, "");
    declareProperty("outputClusterCollection", m_outputClusters, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    //
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
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
      error() << fmt::format("Cannot find energy weighting method {}, choose one from [{}]", m_energyWeight,
                             boost::algorithm::join(weightMethods | boost::adaptors::map_keys, ", "))
              << endmsg;
      return StatusCode::FAILURE;
    }
    weightFunc = it->second;
    // info() << "z_length " << depth << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    //const auto& hits  = *m_inputHits.get();
    const auto& proto = *m_inputProto.get();
    auto& clusters    = *m_outputClusters.createAndPut();

    for (const auto& pcl : proto) {
      auto cl = reconstruct(pcl);

      if (msgLevel(MSG::DEBUG)) {
        debug() << cl.nhits() << " hits: " << cl.energy() / GeV << " GeV, (" << cl.position().x / mm << ", "
                << cl.position().y / mm << ", " << cl.position().z / mm << ")" << endmsg;
      }
      clusters.push_back(cl);
    }

    return StatusCode::SUCCESS;
  }

private:
  eic::Cluster reconstruct(const eic::ConstProtoCluster& pcl) const {
    eic::Cluster cl;
    cl.nhits(pcl.hits_size());

    // no hits
    if (msgLevel(MSG::DEBUG)) {
      debug() << "hit size = " << pcl.hits_size() << endmsg;
    }
    if (pcl.hits_size() == 0) {
      return cl;
    }

    // calculate total energy, find the cell with the maximum energy deposit
    float totalE = 0.;
    float maxE   = 0.;
    // Used to optionally constrain the cluster eta to those of the contributing hits
    float minHitEta = std::numeric_limits<float>::max();
    float maxHitEta = std::numeric_limits<float>::min();
    auto time       = pcl.hits()[0].time();
    auto timeError  = pcl.hits()[0].timeError();
    for (unsigned i = 0; i < pcl.hits().size(); ++i) {
      const auto& hit   = pcl.hits()[i];
      const auto weight = pcl.weights()[i];
      if (msgLevel(MSG::DEBUG)) {
        debug() << "hit energy = " << hit.energy() << " hit weight: " << weight << endmsg;
      }
      auto energy = hit.energy() * weight;
      totalE += energy;
      if (energy > maxE) {
        maxE = energy;
        time = hit.time();
      }
      const float eta = eicd::eta(hit.position());
      if (eta < minHitEta) {
        minHitEta = eta;
      }
      if (eta > maxHitEta) {
        maxHitEta = eta;
      }
    }
    cl.energy(totalE / m_sampFrac);
    cl.energyError(0.);
    cl.time(time);
    cl.timeError(timeError);

    // center of gravity with logarithmic weighting
    float tw = 0.;
    auto v = cl.position();
    for (unsigned i = 0; i < pcl.hits().size(); ++i) {
      const auto& hit   = pcl.hits()[i];
      const auto weight = pcl.weights()[i];
      float w           = weightFunc(hit.energy() * weight, totalE, m_logWeightBase.value(), 0);
      tw += w;
      v = v + (hit.position() * w);
    }
    if (tw == 0.) {
      warning() << "zero total weights encountered, you may want to adjust your weighting parameter." << endmsg;
    }
    cl.position(v / tw);
    cl.positionError({}); // @TODO: Covariance matrix

    // Optionally constrain the cluster to the hit eta values
    if (m_enableEtaBounds) {
      const bool overflow  = (eicd::eta(cl.position()) > maxHitEta);
      const bool underflow = (eicd::eta(cl.position()) < minHitEta);
      if (overflow || underflow) {
        const double newEta           = overflow ? maxHitEta : minHitEta;
        const double newTheta         = eicd::etaToAngle(newEta);
        const double newR             = eicd::magnitude(cl.position());
        const double newPhi           = eicd::angleAzimuthal(cl.position());
        cl.position(eicd::sphericalToVector(newR, newTheta, newPhi));
        if (msgLevel(MSG::DEBUG)) {
          debug() << "Bound cluster position to contributing hits due to " << (overflow ? "overflow" : "underflow")
                  << endmsg;
        }
      }
    }

    // Additional convenience variables

    // best estimate on the cluster direction is the cluster position
    // for simple 2D CoG clustering
    cl.intrinsicTheta(eicd::anglePolar(cl.position()));
    cl.intrinsicPhi(eicd::angleAzimuthal(cl.position()));
    // TODO errors

    // Calculate radius
    // @TODO: add skewness
    if (cl.nhits() > 1) {
      double radius = 0;
      for (const auto& hit : pcl.hits()) {
        const auto delta = cl.position() - hit.position();
        radius += delta * delta;
      }
      radius = sqrt((1. / (cl.nhits() - 1.)) * radius);
      cl.addshapeParameters(radius);
      cl.addshapeParameters(0 /* skewness */); // skewness not yet calculated
    }

    // Optionally store the MC truth associated with the first hit in this cluster
    // FIXME no connection between cluster and truth in edm4hep
    // if (mcHits) {
    //  const auto& mc_hit    = (*mcHits)[pcl.hits(0).ID.value];
    //  cl.mcID({mc_hit.truth().trackID, m_kMonteCarloSource});
    //}

    return cl;
  }
};

DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

