//
//  GeoSvc.cxx
//
//
//  Created by Julia Hrdinka on 30/03/15.
//
//

#include "GeoSvc.h"
#include "GaudiKernel/Service.h"
//#include "GeoConstruction.h"
#include "TGeoManager.h"

#include "DD4hep/Printout.h"

#include "JugBase/ACTSLogger.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
using namespace Gaudi;

DECLARE_COMPONENT(GeoSvc)

GeoSvc::GeoSvc(const std::string& name, ISvcLocator* svc)
    : base_class(name, svc)
    //, m_incidentSvc("IncidentSvc", "GeoSvc")
    , m_trackingGeo(nullptr)
    , m_dd4hepgeo(0)
    //, m_geant4geo(0)
    , m_log(msgSvc(), name) {}

GeoSvc::~GeoSvc() {
  if (m_dd4hepgeo) {
    try {
      m_dd4hepgeo->destroyInstance();
      m_dd4hepgeo = 0;
    } catch (...) {
    }
  }
}

StatusCode GeoSvc::initialize() {
  StatusCode sc = Service::initialize();
  if (!sc.isSuccess())
    return sc;
  // Turn off TGeo printouts if appropriate for the msg level
  if (msgLevel() >= MSG::INFO) {
    TGeoManager::SetVerboseLevel(0);
  }
  uint printoutLevel = msgLevel();
  dd4hep::setPrintLevel(dd4hep::PrintLevel(printoutLevel));
  // m_incidentSvc->addListener(this, "GeometryFailure");
  if (buildDD4HepGeo().isFailure())
    m_log << MSG::ERROR << "Could not build DD4Hep geometry" << endmsg;
  else
    m_log << MSG::INFO << "DD4Hep geometry SUCCESSFULLY built" << endmsg;

  // if (buildGeant4Geo().isFailure())
  //  m_log << MSG::ERROR << "Could not build Geant4 geometry" << endmsg;
  // else
  //  m_log << MSG::INFO << "Geant4 geometry SUCCESSFULLY built" << endmsg;
  // if (m_failureFlag) {
  //  return StatusCode::FAILURE;
  //}
  Acts::Logging::Level geoMsgLevel;
  switch (msgLevel()) {
  case (MSG::DEBUG):
    geoMsgLevel = Acts::Logging::DEBUG;
    break;
  case (MSG::VERBOSE):
    geoMsgLevel = Acts::Logging::VERBOSE;
    break;
  case (MSG::INFO):
    geoMsgLevel = Acts::Logging::INFO;
    break;
  case (MSG::WARNING):
    geoMsgLevel = Acts::Logging::WARNING;
    break;
  case (MSG::FATAL):
    geoMsgLevel = Acts::Logging::FATAL;
    break;
  case (MSG::ERROR):
    geoMsgLevel = Acts::Logging::ERROR;
    break;
  default:
    geoMsgLevel = Acts::Logging::VERBOSE;
  }
  m_trackingGeo = std::move(Acts::convertDD4hepDetector(
      m_dd4hepgeo->world(), geoMsgLevel, Acts::equidistant, Acts::equidistant,
      Acts::equidistant));
  return StatusCode::SUCCESS;
}

StatusCode GeoSvc::finalize() { return StatusCode::SUCCESS; }

StatusCode GeoSvc::buildDD4HepGeo() {
  // we retrieve the the static instance of the DD4HEP::Geometry
  m_dd4hepgeo = &(dd4hep::Detector::getInstance());
  m_dd4hepgeo->addExtension<IGeoSvc>(this);

  // load geometry
  for (auto& filename : m_xmlFileNames) {
    m_log << MSG::INFO << "loading geometry from file:  '" << filename << "'" << endmsg;
    m_dd4hepgeo->fromCompact(filename);
  }
  m_dd4hepgeo->volumeManager();
  m_dd4hepgeo->apply("DD4hepVolumeManager", 0, 0);

  return StatusCode::SUCCESS;
}

dd4hep::Detector* GeoSvc::detector() { return (m_dd4hepgeo); }

dd4hep::DetElement GeoSvc::getDD4HepGeo() { return (detector()->world()); }

//StatusCode GeoSvc::buildGeant4Geo() {
//  std::shared_ptr<G4VUserDetectorConstruction> detector(new det::GeoConstruction(*lcdd()));
//  m_geant4geo = detector;
//  if (m_geant4geo) {
//    return StatusCode::SUCCESS;
//  } else
//    return StatusCode::FAILURE;
//}

//G4VUserDetectorConstruction* GeoSvc::getGeant4Geo() { return (m_geant4geo.get()); }

//void GeoSvc::handle(const Incident& inc) {
//  error() << "Handling incident '" << inc.type() << "'" << endmsg;
//  if (!inc.type().compare("GeometryFailure")) {
//    m_failureFlag = true;
//  }
//}
