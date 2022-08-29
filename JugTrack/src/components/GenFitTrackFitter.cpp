// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten

#include "GenFitTrackFitter.h"
// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "JugBase/BField/DD4hepBField.h"
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugTrack/GeometryContainers.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/Track.hpp"

#include "edm4eic/TrackerHitCollection.h"
#include "edm4eic/vector_utils.h"

#include <functional>
#include <random>
#include <stdexcept>
#include <vector>

//# genfit
#include "ConstField.h"
#include "DAF.h"
#include "Exception.h"
#include "FieldManager.h"
#include "KalmanFitterRefTrack.h"
#include "MaterialEffects.h"
#include "RKTrackRep.h"
#include "StateOnPlane.h"
#include "TGeoMaterialInterface.h"
#include "Track.h"
#include "TrackPoint.h"
//#include <EventDisplay.h>
#include "HelixTrackModel.h"
#include "PlanarMeasurement.h"
//#include "MeasurementCreator.h"
#include "WireMeasurement.h"

#include "TDatabasePDG.h"
#include "TEveManager.h"
#include "TGeoManager.h"
#include "TRandom.h"
#include "TVector3.h"
#include <vector>

namespace Jug::Reco {

GenFitTrackFitter::GenFitTrackFitter(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
  declareProperty("inputHitCollection", m_inputHitCollection, "");
  declareProperty("initialTrackParameters", m_initialTrackParameters, "");
  declareProperty("inputProtoTracks", m_inputProtoTracks, "");
  declareProperty("trackParameters", m_foundTracks, "");
  declareProperty("outputTrajectories", m_outputTrajectories, "");
}

StatusCode GenFitTrackFitter::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  genfit::FieldManager::getInstance()->init(new FieldImp(m_geoSvc->detector()));
  // 0., 0., m_geoSvc->centralMagneticField() * 10.0)); // gentfit uses kilo-Gauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  // copy the whole map to get around genfit's interface
  // this should be returning a const&
  m_detPlaneMap = m_geoSvc->getDetPlaneMap();
  m_surfaceMap  = m_geoSvc->getDD4hepSurfaceMap();

  return StatusCode::SUCCESS;
}

StatusCode GenFitTrackFitter::execute() {
  // Read input data
  const edm4eic::TrackerHitCollection* hits            = m_inputHitCollection.get();
  const TrackParametersContainer* initialParameters = m_initialTrackParameters.get();
  const ProtoTrackContainer* protoTracks            = m_inputProtoTracks.get();

  // TrajectoryContainer trajectories;
  // commented out unused variables
  /*auto trajectories    =*/m_outputTrajectories.createAndPut();
  /*auto trackParameters =*/m_foundTracks.createAndPut();

  int n_tracks       = initialParameters->size();
  int n_proto_tracks = protoTracks->size();
  // Unused variable
  // int ID             = 0;

  // Assuming init track parameters have been match with proto tracks by index
  if (n_proto_tracks != n_tracks) {
    warning() << " Number of proto tracks does not match the initial track parameters." << endmsg;
  }

  for (int itrack = 0; itrack < std::min(n_tracks, n_proto_tracks); itrack++) {
    const auto& track_param = (*initialParameters)[itrack];
    const auto& proto_track = (*protoTracks)[itrack];

    if (msgLevel(MSG::DEBUG)) {
      debug() << "track mom : " << track_param.absoluteMomentum() << endmsg;
    }
    if (hits->size() < 2) {
      return StatusCode::SUCCESS;
    }
    // init fitter
    // genfit::KalmanFitterRefTrack fitter;
    // fitter.setDebugLvl(1);
    genfit::DAF fitter;
    // genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

    ROOT::Math::XYZVector tp(track_param.momentum()[0], track_param.momentum()[1], track_param.momentum()[2]);
    auto first_hit       = (*hits)[proto_track[0]];
    auto first_hit_phi   = edm4eic::angleAzimuthal(first_hit.getPosition());
    auto track_param_phi = tp.phi();
    if (msgLevel(MSG::DEBUG)) {
      debug() << " first hit phi:  " << first_hit_phi << endmsg;
      debug() << "init track phi:  " << track_param_phi << endmsg;
    }
    if (std::fabs(first_hit_phi - track_param_phi) > 0.15) {
      warning() << "Seed directions does not match first hit phi. " << endmsg;
      continue;
    }

    // start values for the fit, e.g. from pattern recognition
    TVector3 pos(0, 0, 0);
    TVector3 mom(track_param.momentum()[0], track_param.momentum()[1], track_param.momentum()[2]);
    TMatrixDSym covM(6);
    covM(0, 0) = 0.001;
    covM(1, 1) = 0.001;
    covM(2, 2) = 1.0;
    covM(3, 3) = 0.05 * track_param.momentum()[0] * 0.05 * track_param.momentum()[0];
    covM(4, 4) = 0.05 * track_param.momentum()[1] * 0.05 * track_param.momentum()[1];
    covM(5, 5) = 0.05 * track_param.momentum()[2] * 0.05 * track_param.momentum()[2];

    if (msgLevel(MSG::DEBUG)) {
      debug() << "covM = " << covM(0, 0) << "," << covM(1, 1) << "," << covM(2, 2) << "," << covM(3, 3) << ","
              << covM(4, 4) << "," << covM(5, 5) << " " << endmsg;
    }

    // trackrep
    // @FIXME: raw new should be avoided, either place on the stack or use
    // std::unique_ptr<>
    genfit::AbsTrackRep* electron_rep = new genfit::RKTrackRep(11);
    // unusud
    // genfit::AbsTrackRep* positron_rep = new genfit::RKTrackRep(-11);
    // genfit::AbsTrackRep* piplus_rep = new genfit::RKTrackRep(211);
    // genfit::AbsTrackRep* piminus_rep = new genfit::RKTrackRep(-211);

    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(electron_rep);
    stateSmeared.setPosMomCov(pos, mom, covM);

    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    stateSmeared.get6DStateCov(seedState, seedCov);
    // genfit::Track fitTrack(rep, seedState, seedCov);

    // create track
    genfit::Track fitTrack(electron_rep, seedState, seedCov);
    // genfit::Track fitTrack(electron_rep, pos, mom);
    // fitTrack.addTrackRep(positron_rep);
    // fitTrack.addTrackRep(piplus_rep );
    // fitTrack.addTrackRep(piminus_rep);

    if (msgLevel(MSG::DEBUG)) {
      debug() << (*hits).size() << " hits " << endmsg;
    }

    int nhit = 0;
    for (int ihit : proto_track) {
      const auto& ahit = (*hits)[ihit];

      const auto* vol_ctx        = m_geoSvc->cellIDPositionConverter()->findContext(ahit.getCellID());
      auto vol_id         = vol_ctx->identifier;
      auto volman         = m_geoSvc->detector()->volumeManager();
      auto alignment      = volman.lookupDetElement(vol_id).nominal();
      auto local_position = alignment.worldToLocal(
          {ahit.getPosition().x / 10.0, ahit.getPosition().y / 10.0, ahit.getPosition().z / 10.0});
      auto* surf = m_surfaceMap[vol_id];
      auto local_position2 =
          surf->globalToLocal({ahit.getPosition().x / 10.0, ahit.getPosition().y / 10.0, ahit.getPosition().z / 10.0});

      TMatrixDSym hitCov(2);
      hitCov.UnitMatrix();
      hitCov(0, 0) = ahit.getPositionError().xx / (100.0); // go from mm^2 to  cm^2
      hitCov(1, 1) = ahit.getPositionError().yy / (100.0); // go from mm^2 to  cm^2

      if (msgLevel(MSG::DEBUG)) {
        debug() << "------------------------------------ " << endmsg;
        debug() << " hit position     : " << ahit.getPosition().x / 10 << " " << ahit.getPosition().y / 10 << " "
                << ahit.getPosition().z / 10 << endmsg;
        debug() << " dd4hep loc  pos  : " << local_position.x() << " " << local_position.y() << " "
                << local_position.z() << endmsg;
        debug() << " dd4hep surf pos  : " << local_position2.u() << " " << local_position2.v() << endmsg;
      }

      /** \todo Add check for XZ segmentations to use the right local coordinates.
       *  Unlike acts, the conversion to the local system isn't going from 3D -> 2D.
       *  Thefore there is one coordinate that is zero. Which one depends on the the segmentation
       *  type XY, XZ, etc. For XY the Z-coordinate is zero. For the XZ, the Y-coordinate is zero.
       */
      TVectorD hitCoords(2);
      hitCoords[0] = local_position2.u();
      hitCoords[1] = local_position2.v();
      if (msgLevel(MSG::DEBUG)) {
        debug() << "covariance matrix :  " << hitCov(0, 0) << " " << hitCov(1, 1) << " " << endmsg;
        debug() << "  hit coordinates :  " << hitCoords[0] << " " << hitCoords[1] << " " << endmsg;
      }
      auto* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 1 /** type **/, nhit, nullptr);

      // measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(point, u_dir, v_dir)),
      measurement->setPlane(m_detPlaneMap[vol_id], vol_id);
      fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
      // positronFitTrack.insertPoint(new genfit::TrackPoint(measurement, &positronFitTrack));

      nhit++;
    }

    // do the fit
    try {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "Electron track: " << endmsg;
        fitTrack.checkConsistency();
      }
      fitter.processTrack(&fitTrack, true);
      bool isConverged = fitTrack.getFitStatus()->isFitConverged();
      if (!isConverged) {
        fitter.processTrack(&fitTrack, true);
      }

      isConverged = fitTrack.getFitStatus()->isFitConverged();
      if (!isConverged) {
        fitter.processTrack(&fitTrack);
      }

      // print fit result
      fitTrack.getFittedState().Print();
      // isConverged = fitTrack.getFitStatus()->isFitConverged();

      bool isFitted = fitTrack.getFitStatus()->isFitted();
      float chi2    = fitTrack.getFitStatus()->getChi2();
      float ndf     = fitTrack.getFitStatus()->getNdf();
      // unused
      // float charge      = fitTrack.getFitStatus()->getCharge();

      TVector3 vertexPos;
      TVector3 vertexMom;
      TMatrixDSym vertexCov;
      genfit::MeasuredStateOnPlane state = fitTrack.getFittedState(); // copy
      TVector3 vertex(0, 0, 0);
      TVector3 axis(0, 0, 1);
      // state.extrapolateToPoint(vertex);
      // or alternatively
      state.extrapolateToLine(vertex, axis);
      state.getPosMomCov(vertexPos, vertexMom, vertexCov);

      if (msgLevel(MSG::DEBUG)) {
        debug() << "Electron track: " << endmsg;
        fitTrack.checkConsistency();
        debug() << "vertex pos: " << vertexPos.x() << ", " << vertexPos.y() << ", " << vertexPos.z() << endmsg;
        debug() << "vertex mom: " << vertexMom.x() << ", " << vertexMom.y() << ", " << vertexMom.z() << endmsg;
        debug() << "track status: " << endmsg;
        debug() << " fitted    = " << isFitted << endmsg;
        debug() << " converged =" << isConverged << endmsg;
        debug() << " chi2/ndf    = " << isFitted << "/" << ndf << " = " << chi2 / ndf << endmsg;
        debug() << " charge =" << isConverged << endmsg;
        // debug() << "Positron track: " << endmsg;
        // positronFitTrack.checkConsistency();
      }
    } catch (genfit::Exception& e) {
      warning() << e.what() << endmsg;
      warning() << "Exception, next track" << endmsg;
      continue;
    }

    // edm4eic::TrackParameters electron_track_params({ID++, algorithmID()}, {0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},

    // TrackParameters(edm4eic::Index ID, edm4eic::FloatPair loc, edm4eic::FloatPair locError, edm4eic::Direction direction,
    // edm4eic::Direction directionError, float qOverP, float qOverPError, float time, float timeError);
    // tracks->push_back(electron_track_params);

    // delete fitter;
  }

  return StatusCode::SUCCESS;
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(GenFitTrackFitter)
} // namespace Jug::Reco
