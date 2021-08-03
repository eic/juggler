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

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugTrack/BField.h"
#include "JugTrack/GeometryContainers.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/Track.hpp"

#include "eicd/TrackerHitCollection.h"

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

  GenFitTrackFitter::GenFitTrackFitter(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
  {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("initialTrackParameters", m_initialTrackParameters, "");
    declareProperty("inputProtoTracks", m_inputProtoTracks, "");
    declareProperty("foundTracks", m_foundTracks, "");
    declareProperty("outputTrajectories", m_outputTrajectories, "");
  }

  StatusCode GenFitTrackFitter::initialize()
  {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }

    // genfit::FieldManager::getInstance()->init(new genfit::ConstField(
    //    0., 0., m_geoSvc->centralMagneticField() * 10.0)); // gentfit uses kilo-Gauss
    // genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

    return StatusCode::SUCCESS;
  }

  StatusCode GenFitTrackFitter::execute()
  {
    // Read input data
    const eic::TrackerHitCollection* hits              = m_inputHitCollection.get();
    const TrackParametersContainer*  initialParameters = m_initialTrackParameters.get();
    const ProtoTrackContainer*       protoTracks       = m_inputProtoTracks.get();

    const auto& single_track_param = (*initialParameters)[0];
    // TrajectoryContainer trajectories;
    // auto trajectories = m_outputTrajectories.createAndPut();

    // copy the map to get around genfit's interface
    // this should be returning a const&
    auto detPlaneMap = m_geoSvc->getDetPlaneMap();

    debug() << "Single track mom : " << single_track_param.absoluteMomentum() << endmsg;
    if (hits->size() < 2) {
      return StatusCode::SUCCESS;
    }
    // init fitter
    genfit::KalmanFitterRefTrack fitter;
    // genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

    // particle pdg code; pion hypothesis
    const int pdg = 11;

    // start values for the fit, e.g. from pattern recognition
    TVector3 pos(0, 0, 0);
    TVector3 mom(single_track_param.momentum()[0], single_track_param.momentum()[1], single_track_param.momentum()[2]);

    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // create track
    genfit::Track fitTrack(rep, pos, mom);

    debug() << (*hits).size() << " hits " << endmsg;
    int ihit = 0;
    for (const auto& ahit : *hits) {

      // auto volman         = m_geoSvc->detector()->volumeManager();
      // auto alignment      = volman.lookupDetElement(vol_id).nominal();
      // auto local_position = alignment.worldToLocal(global_position);
      auto        vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
      auto        vol_id  = vol_ctx->identifier;
      TMatrixDSym hitCov(2);
      hitCov.UnitMatrix();
      hitCov(0, 0) = ahit.covMatrix().covsym_xx * ahit.covMatrix().covsym_xx;
      hitCov(1, 1) = ahit.covMatrix().covsym_yy * ahit.covMatrix().covsym_yy;

      TVector3 point = {ahit.position().x, ahit.position().y, ahit.position().z};
      TVector3 u_dir = {1, 0, 0};
      u_dir.SetPhi(point.Phi());
      TVector3 v_dir = {0, 0, 1};

      // add some planar hits to track with coordinates I just made up
      TVectorD hitCoords(2);
      hitCoords[0]     = 0;
      hitCoords[1]     = 0;
      auto measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, ahit.cellID(), ihit, nullptr);

      // measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(point, u_dir, v_dir)),
      measurement->setPlane(detPlaneMap[vol_id], ihit);
      fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

      ihit++;
    }
    // check
    fitTrack.checkConsistency();

    // do the fit
    fitter.processTrack(&fitTrack);

    // print fit result
    fitTrack.getFittedState().Print();

    // check
    fitTrack.checkConsistency();

    // delete fitter;

    return StatusCode::SUCCESS;
  }

  DECLARE_COMPONENT(GenFitTrackFitter)
} // namespace Jug::Reco
