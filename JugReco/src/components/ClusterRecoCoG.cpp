// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Chao, Chao Peng, Whitney Armstrong

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
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/MCRecoClusterParticleAssociationCollection.h"
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
 *  Logarithmic weighting is used for mimicking energy deposit in transverse direction
 *
 * \ingroup reco
 */
class ClusterRecoCoG : public GaudiAlgorithm {
private:
  Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};
  Gaudi::Property<double> m_logWeightBase{this, "logWeightBase", 3.6};
  Gaudi::Property<double> m_depthCorrection{this, "depthCorrection", 0.0};
  Gaudi::Property<std::string> m_energyWeight{this, "energyWeight", "log"};
  Gaudi::Property<std::string> m_moduleDimZName{this, "moduleDimZName", ""};
  // Constrain the cluster position eta to be within
  // the eta of the contributing hits. This is useful to avoid edge effects
  // for endcaps.
  Gaudi::Property<bool> m_enableEtaBounds{this, "enableEtaBounds", false};

  DataHandle<eicd::ProtoClusterCollection> m_inputProto{"inputProtoClusterCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::ClusterCollection> m_outputClusters{"outputClusterCollection", Gaudi::DataHandle::Writer, this};

  // Collection for MC hits when running on MC
  Gaudi::Property<std::string> m_mcHits{this, "mcHits", ""};
  // Optional handle to MC hits
  std::unique_ptr<DataHandle<edm4hep::SimCalorimeterHitCollection>> m_mcHits_ptr;

  // Collection for associations when running on MC
  Gaudi::Property<std::string> m_outputAssociations{this, "outputAssociations", ""};
  // Optional handle to MC hits
  std::unique_ptr<DataHandle<eicd::MCRecoClusterParticleAssociationCollection>> m_outputAssociations_ptr;

  // Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  double m_depthCorr{0};
  std::function<double(double, double, double, int)> weightFunc;

public:
  ClusterRecoCoG(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputProtoClusterCollection", m_inputProto, "");
    declareProperty("outputClusterCollection", m_outputClusters, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    // Initialize the optional MC input hit collection if requested
    if (m_mcHits != "") {
      m_mcHits_ptr =
        std::make_unique<DataHandle<edm4hep::SimCalorimeterHitCollection>>(m_mcHits, Gaudi::DataHandle::Reader,
        this);
    }

    // Initialize the optional association collection if requested
    if (m_outputAssociations != "") {
      m_outputAssociations_ptr =
        std::make_unique<DataHandle<eicd::MCRecoClusterParticleAssociationCollection>>(m_outputAssociations, Gaudi::DataHandle::Writer,
        this);
    }

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
    const auto& proto  = *m_inputProto.get();
    auto& clusters     = *m_outputClusters.createAndPut();

    // Optional input MC data
    const edm4hep::SimCalorimeterHitCollection* mchits = nullptr;
    if (m_mcHits_ptr) {
      mchits = m_mcHits_ptr->get();
    }

    // Optional output associations
    eicd::MCRecoClusterParticleAssociationCollection* associations = nullptr;
    if (m_outputAssociations_ptr) {
      associations = m_outputAssociations_ptr->createAndPut();
    }

    for (const auto& pcl : proto) {
      auto cl = reconstruct(pcl);

      if (msgLevel(MSG::DEBUG)) {
        debug() << cl.getNhits() << " hits: " << cl.getEnergy() / GeV << " GeV, (" << cl.getPosition().x / mm << ", "
                << cl.getPosition().y / mm << ", " << cl.getPosition().z / mm << ")" << endmsg;
      }
      clusters.push_back(cl);

      // If mcHits are available, associate cluster with MCParticle
      // 1. find proto-cluster hit with largest energy deposition
      // 2. find first mchit with same CellID
      // 3. assign mchit's MCParticle as cluster truth
      if (m_mcHits_ptr.get() != nullptr && m_outputAssociations_ptr.get() != nullptr) {

        // 1. find pclhit with largest energy deposition
        auto pclhits = pcl.getHits();
        auto pclhit = std::max_element(
          pclhits.begin(),
          pclhits.end(),
          [](const auto& pclhit1, const auto& pclhit2) {
            return pclhit1.getEnergy() < pclhit2.getEnergy();
          }
        );

        // 2. find mchit with same CellID
        // find_if not working, https://github.com/AIDASoft/podio/pull/273
        //auto mchit = std::find_if(
        //  mchits.begin(),
        //  mchits.end(),
        //  [&pclhit](const auto& mchit1) {
        //    return mchit1.getCellID() == pclhit->getCellID();
        //  }
        //);
        auto mchit = mchits->begin();
        for ( ; mchit != mchits->end(); ++mchit) {
          // break loop when CellID match found
          if (mchit->getCellID() == pclhit->getCellID()) {
            break;
          }
        }
        if (!(mchit != mchits->end())) {
          // break if no matching hit found for this CellID
          warning() << "Proto-cluster has highest energy in CellID " << pclhit->getCellID()
                    << ", but no mc hit with that CellID was found." << endmsg;
          break;
        }

        // 3. find mchit's MCParticle
        const auto& mcp = mchit->getContributions(0).getParticle();

        // debug output
        if (msgLevel(MSG::DEBUG)) {
          debug() << "cluster has largest energy in cellID: " << pclhit->getCellID() << endmsg;
          debug() << "pcl hit with highest energy " << pclhit->getEnergy() << " at index " << pclhit->getObjectID().index << endmsg;
          debug() << "corresponding mc hit energy " << mchit->getEnergy() << " at index " << mchit->getObjectID().index << endmsg;
          debug() << "from MCParticle index " << mcp.getObjectID().index << ", PDG " << mcp.getPDG() << ", " << eicd::magnitude(mcp.getMomentum()) << endmsg;
        }

        // set association
        eicd::MutableMCRecoClusterParticleAssociation clusterassoc;
        clusterassoc.setSimID(mcp.getObjectID().index);
        clusterassoc.setWeight(1.0);
        clusterassoc.setRec(cl);
        associations->push_back(clusterassoc);
      } else {
        if (msgLevel(MSG::DEBUG)) {
          debug() << "No mcHitCollection was provided, so no truth association will be performed." << endmsg;
        }
      }
    }

    return StatusCode::SUCCESS;
  }

private:
  eicd::MutableCluster reconstruct(const eicd::ProtoCluster& pcl) const {
    eicd::MutableCluster cl;
    cl.setNhits(pcl.hits_size());

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
    auto time       = pcl.getHits()[0].getTime();
    auto timeError  = pcl.getHits()[0].getTimeError();
    for (unsigned i = 0; i < pcl.getHits().size(); ++i) {
      const auto& hit   = pcl.getHits()[i];
      const auto weight = pcl.getWeights()[i];
      if (msgLevel(MSG::DEBUG)) {
        debug() << "hit energy = " << hit.getEnergy() << " hit weight: " << weight << endmsg;
      }
      auto energy = hit.getEnergy() * weight;
      totalE += energy;
      if (energy > maxE) {
      }
      const float eta = eicd::eta(hit.getPosition());
      if (eta < minHitEta) {
        minHitEta = eta;
      }
      if (eta > maxHitEta) {
        maxHitEta = eta;
      }
    }
    cl.setEnergy(totalE / m_sampFrac);
    cl.setEnergyError(0.);
    cl.setTime(time);
    cl.setTimeError(timeError);

    // center of gravity with logarithmic weighting
    float tw = 0.;
    auto v   = cl.getPosition();
    for (unsigned i = 0; i < pcl.getHits().size(); ++i) {
      const auto& hit   = pcl.getHits()[i];
      const auto weight = pcl.getWeights()[i];
      float w           = weightFunc(hit.getEnergy() * weight, totalE, m_logWeightBase.value(), 0);
      tw += w;
      v = v + (hit.getPosition() * w);
    }
    if (tw == 0.) {
      warning() << "zero total weights encountered, you may want to adjust your weighting parameter." << endmsg;
    }
    cl.setPosition(v / tw);
    cl.setPositionError({}); // @TODO: Covariance matrix

    // Optionally constrain the cluster to the hit eta values
    if (m_enableEtaBounds) {
      const bool overflow  = (eicd::eta(cl.getPosition()) > maxHitEta);
      const bool underflow = (eicd::eta(cl.getPosition()) < minHitEta);
      if (overflow || underflow) {
        const double newEta   = overflow ? maxHitEta : minHitEta;
        const double newTheta = eicd::etaToAngle(newEta);
        const double newR     = eicd::magnitude(cl.getPosition());
        const double newPhi   = eicd::angleAzimuthal(cl.getPosition());
        cl.setPosition(eicd::sphericalToVector(newR, newTheta, newPhi));
        if (msgLevel(MSG::DEBUG)) {
          debug() << "Bound cluster position to contributing hits due to " << (overflow ? "overflow" : "underflow")
                  << endmsg;
        }
      }
    }

    // Additional convenience variables

    // best estimate on the cluster direction is the cluster position
    // for simple 2D CoG clustering
    cl.setIntrinsicTheta(eicd::anglePolar(cl.getPosition()));
    cl.setIntrinsicPhi(eicd::angleAzimuthal(cl.getPosition()));
    // TODO errors

    // Calculate radius
    // @TODO: add skewness
    if (cl.getNhits() > 1) {
      double radius = 0;
      for (const auto& hit : pcl.getHits()) {
        const auto delta = cl.getPosition() - hit.getPosition();
        radius += delta * delta;
      }
      radius = sqrt((1. / (cl.getNhits() - 1.)) * radius);
      cl.addToShapeParameters(radius);
      cl.addToShapeParameters(0 /* skewness */); // skewness not yet calculated
    }

    return cl;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco
