// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Alex Jentsch, Sylvester Joosten, Wouter Deconinck

#include <algorithm>
#include <cmath>
#include <fmt/format.h>

#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"

// Event Model related classes
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/TrackerHitCollection.h"
#include <eicd/vector_utils.h>

namespace Jug::Reco {

class FarForwardParticlesOMD : public GaudiAlgorithm {
private:
  DataHandle<eicd::TrackerHitCollection> m_inputHitCollection{"FarForwardTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::ReconstructedParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer,
                                                                     this};

  //----- Define constants here ------

  Gaudi::Property<double> local_x_offset_station_1{this, "localXOffsetSta1", -762.5006104};
  Gaudi::Property<double> local_x_offset_station_2{this, "localXOffsetSta2", -881.9621277};
  Gaudi::Property<double> local_x_slope_offset{this, "localXSlopeOffset", -59.73075865};
  Gaudi::Property<double> local_y_slope_offset{this, "localYSlopeOffset", 0.0012755};
  Gaudi::Property<double> crossingAngle{this, "crossingAngle", -0.025};
  Gaudi::Property<double> nomMomentum{this, "beamMomentum", 137.5}; // This number is set to 50% maximum beam momentum

  const double aXOMD[2][2] = {{1.6229248, 12.9519653}, {-2.86056525, 0.1830292}};
  const double aYOMD[2][2] = {{0.0000185, -28.599739}, {0.00000925, -2.8795791}};

  double aXOMDinv[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double aYOMDinv[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

public:
  FarForwardParticlesOMD(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputCollection", m_inputHitCollection, "FarForwardTrackerHits");
    declareProperty("outputCollection", m_outputParticles, "ReconstructedParticles");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    double det = aXOMD[0][0] * aXOMD[1][1] - aXOMD[0][1] * aXOMD[1][0];

    if (det == 0) {
      error() << "Reco matrix determinant = 0!"
              << "Matrix cannot be inverted! Double-check matrix!" << endmsg;
      return StatusCode::FAILURE;
    }

    aXOMDinv[0][0] = aXOMD[1][1] / det;
    aXOMDinv[0][1] = -aXOMD[0][1] / det;
    aXOMDinv[1][0] = -aXOMD[1][0] / det;
    aXOMDinv[1][1] = aXOMD[0][0] / det;

    det            = aYOMD[0][0] * aYOMD[1][1] - aYOMD[0][1] * aYOMD[1][0];
    aYOMDinv[0][0] = aYOMD[1][1] / det;
    aYOMDinv[0][1] = -aYOMD[0][1] / det;
    aYOMDinv[1][0] = -aYOMD[1][0] / det;
    aYOMDinv[1][1] = aYOMD[0][0] / det;

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    const eicd::TrackerHitCollection* rawhits = m_inputHitCollection.get();
    auto& rc                                 = *(m_outputParticles.createAndPut());

    // for (const auto& part : mc) {
    //    if (part.genStatus() > 1) {
    //        if (msgLevel(MSG::DEBUG)) {
    //            debug() << "ignoring particle with genStatus = " << part.genStatus() << endmsg;
    //        }
    //        continue;
    //    }

    //---- begin Roman Pot Reconstruction code ----

    int eventReset = 0; // counter for IDing at least one hit per layer
    std::vector<double> hitx;
    std::vector<double> hity;
    std::vector<double> hitz;

    for (const auto& h : *rawhits) {

      // The actual hit position:
      auto pos0 = h.getPosition();
      // auto mom0 = h.momentum;
      // auto pidCode = h.g4ID;
      auto eDep = h.getEdep();

      if (eDep < 0.00001) {
        continue;
      }

      if (eventReset < 2) {
        hitx.push_back(pos0.x - local_x_offset_station_2);
      } // use station 2 for both offsets since it is used for the reference orbit
      else {
        hitx.push_back(pos0.x - local_x_offset_station_2);
      }

      hity.push_back(pos0.y);
      hitz.push_back(pos0.z);

      eventReset++;
    }

    // NB:
    // This is a "dumb" algorithm - I am just checking the basic thing works with a simple single-proton test.
    // Will need to update and modify for generic # of hits for more complicated final-states.

    if (eventReset == 4) {

      // extract hit, subtract orbit offset – this is to get the hits in the coordinate system of the orbit trajectory
      double XL[2] = {hitx[0], hitx[2]};
      double YL[2] = {hity[0], hity[2]};

      double base = hitz[2] - hitz[0];

      if (base == 0) {
        warning() << "Detector separation = 0!"
                << "Cannot calculate slope!" << endmsg;
        return StatusCode::SUCCESS;
      }

      double Xomd[2] = {XL[1], (1000 * (XL[1] - XL[0]) / (base)) - local_x_slope_offset};
      double Xip[2]  = {0.0, 0.0};
      double Yomd[2] = {YL[1], (1000 * (YL[1] - YL[0]) / (base)) - local_y_slope_offset};
      double Yip[2]  = {0.0, 0.0};

      // use the hit information and calculated slope at the RP + the transfer matrix inverse to calculate the Polar
      // Angle and deltaP at the IP

      for (unsigned i0 = 0; i0 < 2; i0++) {
        for (unsigned i1 = 0; i1 < 2; i1++) {
          Xip[i0] += aXOMDinv[i0][i1] * Xomd[i1];
          Yip[i0] += aYOMDinv[i0][i1] * Yomd[i1];
        }
      }

      // convert polar angles to radians
      double rsx = Xip[1] / 1000.;
      double rsy = Yip[1] / 1000.;

      // calculate momentum magnitude from measured deltaP – using thin lens optics.
      double p    = nomMomentum * (1 + 0.01 * Xip[0]);
      double norm = std::sqrt(1.0 + rsx * rsx + rsy * rsy);

      const float prec[3] = {static_cast<float>(p * rsx / norm), static_cast<float>(p * rsy / norm),
                             static_cast<float>(p / norm)};

      //----- end RP reconstruction code ------

      eicd::MutableReconstructedParticle rpTrack;
      rpTrack.setType(0);
      rpTrack.setMomentum({prec});
      rpTrack.setEnergy(std::hypot(eicd::magnitude(rpTrack.getMomentum()), .938272));
      rpTrack.setReferencePoint({0, 0, 0});
      rpTrack.setCharge(1);
      rpTrack.setMass(.938272);
      rpTrack.setGoodnessOfPID(1.);
      rpTrack.setPDG(2122);
      //rpTrack.covMatrix(); // @TODO: Errors
      rc->push_back(rpTrack);

    } // end enough hits for matrix reco

    return StatusCode::SUCCESS;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(FarForwardParticlesOMD)

} // namespace Jug::Reco
