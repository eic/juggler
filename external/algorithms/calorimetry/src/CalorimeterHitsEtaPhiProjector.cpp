// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Jihee Kim, Sylvester Joosten, Chao Peng, Whitney Armstrong, Wouter Deconinck, Chao Peng

/*
 *  A hits converter to prepare dataset for machine learning
 *  It converts hits with (x, y, z, E) to (E, eta, phi) layer by layer
 *  With a defined grid size and ranges, it merge the hits within one grid and drop-off hits
 * out-of-range The capacity of each layer is fixed (padding with zeros), and the hits with least
 * energies that exceed the capacity will be discarded.
 *
 *  Author: Chao Peng (ANL), Jihee Kim, 07/16/2021
 *
 */
#include <algorithm>
#include <bitset>
#include <unordered_map>

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
#include "JugBase/Utilities/Utils.hpp"

// Event Model related classes
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/vector_utils.h"

using namespace Gaudi::Units;
using Point3D = ROOT::Math::XYZPoint;

struct pair_hash {
  template <class T1, class T2> std::size_t operator()(const std::pair<T1, T2>& pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

namespace Jug::Reco {

/** Calorimeter eta-phi projector
 *
 *  A hits converter to prepare dataset for machine learning
 *  It converts hits with (x, y, z, E) to (E, eta, phi) layer by layer
 *  With a defined grid size and ranges, it merge the hits within one grid and drop-off hits
 * out-of-range The capacity of each layer is fixed (padding with zeros), and the hits with least
 * energies that exceed the capacity will be discarded.
 *
 *
 * \ingroup reco
 */
class CalorimeterHitsEtaPhiProjector : public GaudiAlgorithm {
private:
  Gaudi::Property<std::vector<double>> u_gridSizes{this, "gridSizes", {0.001, 0.001 * rad}};
  DataHandle<edm4eic::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                  this};
  DataHandle<edm4eic::CalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                   this};

  double gridSizes[2]{0.0, 0.0};

public:
  CalorimeterHitsEtaPhiProjector(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputHitCollection", m_outputHitCollection, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    if (u_gridSizes.size() != 2) {
      error() << "Expected 2 values for gridSizes, received " << u_gridSizes.size() << endmsg;
      return StatusCode::FAILURE;
    }
    gridSizes[0] = u_gridSizes.value()[0];
    gridSizes[1] = u_gridSizes.value()[1] / rad;

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // Create output collections
    auto& mhits = *m_outputHitCollection.createAndPut();

    // container
    std::unordered_map<std::pair<int64_t, int64_t>, std::vector<edm4eic::CalorimeterHit>, pair_hash> merged_hits;

    for (const auto h : *m_inputHitCollection.get()) {
      auto bins =
          std::make_pair(static_cast<int64_t>(pos2bin(edm4eic::eta(h.getPosition()), gridSizes[0], 0.)),
                         static_cast<int64_t>(pos2bin(edm4eic::angleAzimuthal(h.getPosition()), gridSizes[1], 0.)));
      merged_hits[bins].push_back(h);
    }

    for (const auto& [bins, hits] : merged_hits) {
      const auto ref = hits.front();
      edm4eic::MutableCalorimeterHit hit;
      hit.setCellID(ref.getCellID());
      // TODO, we can do timing cut to reject noises
      hit.setTime(ref.getTime());
      double r   = edm4eic::magnitude(ref.getPosition());
      double eta = bin2pos(bins.first, gridSizes[0], 0.);
      double phi = bin2pos(bins.second, gridSizes[1], 1.);
      hit.setPosition(edm4eic::sphericalToVector(r, edm4eic::etaToAngle(eta), phi));
      hit.setDimension({static_cast<float>(gridSizes[0]), static_cast<float>(gridSizes[1]), 0.});
      // merge energy
      hit.setEnergy(0.);
      for (const auto& h : hits) {
        hit.setEnergy(hit.getEnergy() + h.getEnergy());
      }
      mhits.push_back(hit);
    }

    return StatusCode::SUCCESS;
  }

  static int64_t pos2bin(double val, double cell, double offset = 0.) {
    return int64_t(std::floor((val + 0.5 * cell - offset) / cell));
  }

  static double bin2pos(int64_t bin, double cell, double offset = 0.) { return bin * cell + offset; }

}; // class CalorimeterHitsEtaPhiProjector

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(CalorimeterHitsEtaPhiProjector)

} // namespace Jug::Reco
