// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Chao Peng, Whitney Armstrong

/*
 *  Reconstruct the cluster with Center of Gravity method
 *  Logarithmic weighting is used for mimicing energy deposit in transverse direction
 *
 *  Author: Chao Peng (ANL), 09/27/2020
 */

#include <algorithms/calorimetry/ClusterRecoCoG.h>

#include <algorithm>
#include <functional>
#include <map>

#include <fmt/format.h>
#include <fmt/ranges.h>

// Event Model related classes
#include "edm4hep/utils/vector_utils.h"

namespace algorithms::calorimetry {
namespace {

  // weighting functions (with place holders for hit energy, total energy, one parameter
  double constWeight(double /*E*/, double /*tE*/, double /*p*/) { return 1.0; }
  double linearWeight(double E, double /*tE*/, double /*p*/) { return E; }
  double logWeight(double E, double tE, double base) {
    return std::max(0., base + std::log(E / tE));
  }

  const std::map<std::string, ClusterRecoCoG::WeightFunc> weightMethods{
      {"none", constWeight},
      {"linear", linearWeight},
      {"log", logWeight},
  };
} // namespace

void ClusterRecoCoG::init() {
  // select weighting method
  std::string ew = m_energyWeight;
  // make it case-insensitive
  std::transform(ew.begin(), ew.end(), ew.begin(), [](char s) { return std::tolower(s); });
  if (!weightMethods.count(ew)) {
    std::vector<std::string> keys;
    std::transform(weightMethods.begin(), weightMethods.end(), std::back_inserter(keys),
                   [](const auto& keyvalue) { return keyvalue.first; });
    raise(fmt::format("Cannot find energy weighting method {}, choose one from {}", m_energyWeight,
                      keys));
  }
  m_weightFunc = weightMethods.at(ew);
  info() << fmt::format("Energy weight method set to: {}", ew) << endmsg;
}

void ClusterRecoCoG::process(const ClusterRecoCoG::Input& input,
                             const ClusterRecoCoG::Output& output) const {
  const auto [proto, opt_simhits] = input;
  auto [clusters, opt_assoc]      = output;

  for (const auto& pcl : *proto) {
    auto cl = reconstruct(pcl);

    if (aboveDebugThreshold()) {
      debug() << cl.getNhits() << " hits: " << cl.getEnergy() / dd4hep::GeV << " GeV, ("
              << cl.getPosition().x / dd4hep::mm << ", " << cl.getPosition().y / dd4hep::mm << ", "
              << cl.getPosition().z / dd4hep::mm << ")" << endmsg;
    }
    clusters->push_back(cl);

    // If mcHits are available, associate cluster with MCParticle
    // 1. find proto-cluster hit with largest energy deposition
    // 2. find first mchit with same CellID
    // 3. assign mchit's MCParticle as cluster truth
    if (opt_simhits && opt_assoc) {

      // 1. find pclhit with largest energy deposition
      auto pclhits = pcl.getHits();
      auto pclhit  = std::max_element(pclhits.begin(), pclhits.end(),
                                     [](const auto& pclhit1, const auto& pclhit2) {
                                       return pclhit1.getEnergy() < pclhit2.getEnergy();
                                     });

      // 2. find mchit with same CellID
      // find_if not working, https://github.com/AIDASoft/podio/pull/273
      // auto mchit = std::find_if(
      //  opt_simhits->begin(),
      //  opt_simhits->end(),
      //  [&pclhit](const auto& mchit1) {
      //    return mchit1.getCellID() == pclhit->getCellID();
      //  }
      //);
      auto mchit = opt_simhits->begin();
      for (; mchit != opt_simhits->end(); ++mchit) {
        // break loop when CellID match found
        if (mchit->getCellID() == pclhit->getCellID()) {
          break;
        }
      }
      if (!(mchit != opt_simhits->end())) {
        // error condition should not happen
        // break if no matching hit found for this CellID
        warning() << "Proto-cluster has highest energy in CellID " << pclhit->getCellID()
                  << ", but no mc hit with that CellID was found." << endmsg;
        info() << "Proto-cluster hits: " << endmsg;
        for (const auto& pclhit1 : pclhits) {
          info() << pclhit1.getCellID() << ": " << pclhit1.getEnergy() << endmsg;
        }
        info() << "MC hits: " << endmsg;
        for (const auto& mchit1 : *opt_simhits) {
          info() << mchit1.getCellID() << ": " << mchit1.getEnergy() << endmsg;
        }
        break;
      }

      // 3. find mchit's MCParticle
      const auto& mcp = mchit->getContributions(0).getParticle();

      // debug output
      if (aboveDebugThreshold()) {
        debug() << "cluster has largest energy in cellID: " << pclhit->getCellID() << endmsg;
        debug() << "pcl hit with highest energy " << pclhit->getEnergy() << " at index "
                << pclhit->getObjectID().index << endmsg;
        debug() << "corresponding mc hit energy " << mchit->getEnergy() << " at index "
                << mchit->getObjectID().index << endmsg;
        debug() << "from MCParticle index " << mcp.getObjectID().index << ", PDG " << mcp.getPDG()
                << ", " << edm4hep::utils::magnitude(mcp.getMomentum()) << endmsg;
      }

      // set association
      edm4eic::MutableMCRecoClusterParticleAssociation clusterassoc;
      clusterassoc.setRecID(cl.getObjectID().index);
      clusterassoc.setSimID(mcp.getObjectID().index);
      clusterassoc.setWeight(1.0);
      clusterassoc.setRec(cl);
      // clusterassoc.setSim(mcp);
      opt_assoc->push_back(clusterassoc);
    } else {
      if (aboveDebugThreshold()) {
        debug() << "No mcHitCollection was provided, so no truth association will be performed."
                << endmsg;
      }
    }
  }
}

edm4eic::MutableCluster ClusterRecoCoG::reconstruct(const edm4eic::ProtoCluster& pcl) const {
  edm4eic::MutableCluster cl;
  cl.setNhits(pcl.hits_size());

  // no hits
  if (aboveDebugThreshold()) {
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
    if (aboveDebugThreshold()) {
      debug() << "hit energy = " << hit.getEnergy() << " hit weight: " << weight << endmsg;
    }
    auto energy = hit.getEnergy() * weight;
    totalE += energy;
    if (energy > maxE) {
    }
    const float eta = edm4hep::utils::eta(hit.getPosition());
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
    float w           = m_weightFunc(hit.getEnergy() * weight, totalE, m_logWeightBase.value());
    tw += w;
    v = v + (hit.getPosition() * w);
  }
  if (tw == 0.) {
    warning() << "zero total weights encountered, you may want to adjust your weighting parameter."
              << endmsg;
  }
  cl.setPosition(v / tw);
  cl.setPositionError({}); // @TODO: Covariance matrix

  // Optionally constrain the cluster to the hit eta values
  if (m_enableEtaBounds) {
    const bool overflow  = (edm4hep::utils::eta(cl.getPosition()) > maxHitEta);
    const bool underflow = (edm4hep::utils::eta(cl.getPosition()) < minHitEta);
    if (overflow || underflow) {
      const double newEta   = overflow ? maxHitEta : minHitEta;
      const double newTheta = edm4hep::utils::etaToAngle(newEta);
      const double newR     = edm4hep::utils::magnitude(cl.getPosition());
      const double newPhi   = edm4hep::utils::angleAzimuthal(cl.getPosition());
      cl.setPosition(edm4hep::utils::sphericalToVector(newR, newTheta, newPhi));
      if (aboveDebugThreshold()) {
        debug() << "Bound cluster position to contributing hits due to "
                << (overflow ? "overflow" : "underflow") << endmsg;
      }
    }
  }

  // Additional convenience variables

  // best estimate on the cluster direction is the cluster position
  // for simple 2D CoG clustering
  cl.setIntrinsicTheta(edm4hep::utils::anglePolar(cl.getPosition()));
  cl.setIntrinsicPhi(edm4hep::utils::angleAzimuthal(cl.getPosition()));
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

} // namespace algorithms::calorimetry
