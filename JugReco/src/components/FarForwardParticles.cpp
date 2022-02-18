#include <algorithm>
#include <cmath>
#include <fmt/format.h>

#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

class FarForwardParticles : public GaudiAlgorithm, AlgorithmIDMixin<> {
  DataHandle<eic::TrackerHitCollection> m_inputHitCollection{"FarForwardTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer,
                                                                     this};

  //----- Define constants here ------

  Gaudi::Property<double> local_x_offset_station_1{this, "localXOffsetSta1", -833.3878326};
  Gaudi::Property<double> local_x_offset_station_2{this, "localXOffsetSta2", -924.342804};
  Gaudi::Property<double> local_x_slope_offset{this, "localXSlopeOffset", -0.00622147};
  Gaudi::Property<double> local_y_slope_offset{this, "localYSlopeOffset", -0.0451035};
  Gaudi::Property<double> crossingAngle{this, "crossingAngle", -0.025};
  Gaudi::Property<double> nomMomentum{this, "beamMomentum", 275.0};

  Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
  Gaudi::Property<std::string> m_readout{this, "readoutClass", ""};
  Gaudi::Property<std::string> m_layerField{this, "layerField", ""};
  Gaudi::Property<std::string> m_sectorField{this, "sectorField", ""};
  SmartIF<IGeoSvc>             m_geoSvc;
  dd4hep::BitFieldCoder*       id_dec = nullptr;
  size_t                       sector_idx, layer_idx;

  Gaudi::Property<std::string>              m_localDetElement{this, "localDetElement", ""};
  Gaudi::Property<std::vector<std::string>> u_localDetFields{this, "localDetFields", {}};
  dd4hep::DetElement                        local;
  size_t                                    local_mask = ~0;



  const double aXRP[2][2] = {{2.102403743, 29.11067626}, {0.186640381, 0.192604619}};
  const double aYRP[2][2] = {{0.0000159900, 3.94082098}, {0.0000079946, -0.1402995}};

  double aXRPinv[2][2] = {0.0, 0.0};
  double aYRPinv[2][2] = {0.0, 0.0};

public:
  FarForwardParticles(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputCollection", m_inputHitCollection, "FarForwardTrackerHits");
    declareProperty("outputCollection", m_outputParticles, "ReconstructedParticles");
  }

  // See Wouter's example to extract local coordinates CalorimeterHitReco.cpp
  //includes DDRec/CellIDPositionConverter.here
  //See tutorial
  // auto converter = m_GeoSvc ....
  //https://eicweb.phy.anl.gov/EIC/juggler/-/blob/master/JugReco/src/components/CalorimeterHitReco.cpp - line 200
  //include the Eigen libraries, used in ACTS, for the linear algebra.

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()){
      return StatusCode::FAILURE;
    }
    m_geoSvc = service(m_geoSvcName);
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
              << endmsg;
      return StatusCode::FAILURE;
    }

    // do not get the layer/sector ID if no readout class provided
    if (m_readout.value().empty()) {
      return StatusCode::SUCCESS;
    }

    auto id_spec = m_geoSvc->detector()->readout(m_readout).idSpec();
    try {
      id_dec = id_spec.decoder();
      if (m_sectorField.value().size()) {
        sector_idx = id_dec->index(m_sectorField);
        info() << "Find sector field " << m_sectorField.value() << ", index = " << sector_idx << endmsg;
      }
      if (m_layerField.value().size()) {
        layer_idx = id_dec->index(m_layerField);
        info() << "Find layer field " << m_layerField.value() << ", index = " << sector_idx << endmsg;
      }
    } catch (...) {
      error() << "Failed to load ID decoder for " << m_readout << endmsg;
      return StatusCode::FAILURE;
    }

    // local detector name has higher priority
    if (m_localDetElement.value().size()) {
      try {
        local = m_geoSvc->detector()->detector(m_localDetElement.value());
        info() << "Local coordinate system from DetElement " << m_localDetElement.value()
                << endmsg;
      } catch (...) {
        error() << "Failed to locate local coordinate system from DetElement "
                << m_localDetElement.value() << endmsg;
        return StatusCode::FAILURE;
      }
    // or get from fields
    } else {
      std::vector<std::pair<std::string, int>> fields;
      for (auto& f : u_localDetFields.value()) {
        fields.push_back({f, 0});
      }
      local_mask = id_spec.get_mask(fields);
      // use all fields if nothing provided
      if (fields.empty()) {
        local_mask = ~0;
      }
      //info() << fmt::format("Local DetElement mask {:#064b} from fields [{}]", local_mask,
      //                      fmt::join(fields, ", "))
      //        << endmsg;
    }
      
    double det = aXRP[0][0] * aXRP[1][1] - aXRP[0][1] * aXRP[1][0];

    if (det == 0) {
      error() << "Reco matrix determinant = 0!"
              << "Matrix cannot be inverted! Double-check matrix!" << endmsg;
      return StatusCode::FAILURE;
    }

    aXRPinv[0][0] = aXRP[1][1] / det;
    aXRPinv[0][1] = -aXRP[0][1] / det;
    aXRPinv[1][0] = -aXRP[1][0] / det;
    aXRPinv[1][1] = aXRP[0][0] / det;

    det           = aYRP[0][0] * aYRP[1][1] - aYRP[0][1] * aYRP[1][0];
    aYRPinv[0][0] = aYRP[1][1] / det;
    aYRPinv[0][1] = -aYRP[0][1] / det;
    aYRPinv[1][0] = -aYRP[1][0] / det;
    aYRPinv[1][1] = aYRP[0][0] / det;

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    const eic::TrackerHitCollection* rawhits = m_inputHitCollection.get();
    auto& rc                                 = *(m_outputParticles.createAndPut());

    auto converter = m_geoSvc->cellIDPositionConverter();

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

    int32_t idx = 0;
    for (const auto& h : *rawhits) {

      auto  cellID   = h.cellID();
      // The actual hit position in Global Coordinates
      //auto pos0 = h.position();
  
      auto gpos = converter->position(cellID);
        // local positions
        if (m_localDetElement.value().empty()) {
          auto volman = m_geoSvc->detector()->volumeManager();
          local       = volman.lookupDetElement(cellID);
        }
        auto pos0 = local.nominal().worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z())); //hit position in local coordinates

  
      // auto mom0 = h.momentum;
      // auto pidCode = h.g4ID;
      auto eDep = h.edep();



      if (eDep < 0.00001) {
        continue;
      }

      if (eventReset < 2) {
        hitx.push_back(pos0.x());// - local_x_offset_station_2);
      } // use station 2 for both offsets since it is used for the reference orbit
      else {
        hitx.push_back(pos0.x()); // - local_x_offset_station_2);
      }

      hity.push_back(pos0.y());
      hitz.push_back(pos0.z());

      eventReset++;
    }

    // NB:
    // This is a "dumb" algorithm - I am just checking the basic thing works with a simple single-proton test.
    // Will need to update and modify for generic # of hits for more complicated final-states.

    if (eventReset == 4) {

      // extract hit, subtract orbit offset – this is to get the hits in the coordinate system of the orbit
      // trajectory
      double XL[2] = {hitx[0], hitx[2]};
      double YL[2] = {hity[0], hity[2]};

      double base = hitz[2] - hitz[0];

      if (base == 0) {
        warning() << "Detector separation = 0!"
                << "Cannot calculate slope!" << endmsg;
        return StatusCode::SUCCESS;
      }

      double Xrp[2], Xip[2] = {0.0, 0.0};
      Xrp[0] = XL[1];
      Xrp[1] = (1000 * (XL[1] - XL[0]) / (base)) - local_x_slope_offset; //- _SX0RP_;
      double Yrp[2], Yip[2] = {0.0, 0.0};
      Yrp[0] = YL[1];
      Yrp[1] = (1000 * (YL[1] - YL[0]) / (base)) - local_y_slope_offset; //- _SY0RP_;

      // use the hit information and calculated slope at the RP + the transfer matrix inverse to calculate the
      // Polar Angle and deltaP at the IP

      for (unsigned i0 = 0; i0 < 2; i0++) {
        for (unsigned i1 = 0; i1 < 2; i1++) {
          Xip[i0] += aXRPinv[i0][i1] * Xrp[i1];
          Yip[i0] += aYRPinv[i0][i1] * Yrp[i1];
        }
      }

      // convert polar angles to radians
      double rsx = Xip[1] / 1000., rsy = Yip[1] / 1000.;

      // calculate momentum magnitude from measured deltaP – using thin lens optics.
      double p    = nomMomentum * (1 + 0.01 * Xip[0]);
      double norm = std::sqrt(1.0 + rsx * rsx + rsy * rsy);

      double prec[3];
      prec[0] = p * rsx / norm;
      prec[1] = p * rsy / norm;
      prec[2] = p / norm;

      // double pT_reco = std::sqrt(prec[0]*prec[0] + prec[1]*prec[1]);
      float p_reco = std::sqrt(prec[0] * prec[0] + prec[1] * prec[1] + prec[2] * prec[2]);

      //----- end RP reconstruction code ------

      eic::VectorXYZ recoMom(prec[0], prec[1], prec[2]);
      eic::VectorXYZ primVtx(0, 0, 0);
      eic::Direction thetaPhi(recoMom.theta(), recoMom.phi());
      eic::Weight wgt(1);

      // ReconstructedParticle (eic::Index ID, eic::VectorXYZ p, eic::VectorXYZ v, float time, std::int32_t pid,
      // std::int16_t status, std::int16_t charge, eic::Weight weight, eic::Direction direction, float momentum,
      // float energy, float mass)

      // eic::ReconstructedParticle rec_part{{part.ID(), algorithmID()},
      // mom3s,
      //{part.vs().x, part.vs().y, part.vs().z},
      // static_cast<float>(part.time()),
      // part.pdgID(),
      // 0,
      // static_cast<int16_t>(part.charge()),
      // 1.,
      //{mom3s.theta(), mom3s.phi()},
      // static_cast<float>(moms),
      // static_cast<float>(Es),
      // static_cast<float>(part.mass())};

      eic::ReconstructedParticle rpTrack;
      rpTrack.ID({idx++, algorithmID()});
      rpTrack.p(recoMom);
      rpTrack.v(primVtx);
      rpTrack.time(0);
      rpTrack.pid(2122);
      rpTrack.status(0);
      rpTrack.charge(1);
      rpTrack.weight(1.);
      rpTrack.direction({recoMom.theta(), recoMom.phi()});
      rpTrack.momentum(p_reco);
      rpTrack.energy(std::hypot(p_reco, .938272));
      rpTrack.mass(.938272);
      rc->push_back(rpTrack);

    } // end enough hits for matrix reco

    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(FarForwardParticles)

} // namespace Jug::Reco

