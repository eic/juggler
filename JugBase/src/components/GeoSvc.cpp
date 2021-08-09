#include "GeoSvc.h"
#include "GaudiKernel/Service.h"
//#include "GeoConstruction.h"
#include "TGeoManager.h"

#include "DD4hep/Printout.h"

#include "JugBase/ACTSLogger.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

// genfit
#include "ConstField.h"
#include "DAF.h"
#include "Exception.h"
#include "FieldManager.h"
#include "KalmanFitterRefTrack.h"
#include "StateOnPlane.h"
#include "Track.h"
#include "TrackPoint.h"
#include "MaterialEffects.h"
#include "RKTrackRep.h"
#include "TGeoMaterialInterface.h"
#include "PlanarMeasurement.h"

void draw_surfaces(std::shared_ptr<const Acts::TrackingGeometry> trk_geo, const std::string& fname)
{
  using namespace Acts;
  Acts::GeometryContext tgContext = Acts::GeometryContext();
  std::vector<const Surface*> surfaces;

  trk_geo->visitSurfaces([&](const Acts::Surface* surface) {
    // for now we just require a valid surface
    if (not surface) {
      std::cout << " Not a surface \n";
      return;
    }
    surfaces.push_back(surface);
  });
  std::ofstream os;
  os.open(fname);
  os << std::fixed << std::setprecision(6);
  size_t nVtx = 0;
  for (const auto& srfx : surfaces) {
    const PlaneSurface*                 srf    = dynamic_cast<const PlaneSurface*>(srfx);
    const PlanarBounds*                 bounds = dynamic_cast<const PlanarBounds*>(&srf->bounds());
    for (const auto& vtxloc : bounds->vertices()) {
      Vector3 vtx = srf->transform(tgContext) * Vector3(vtxloc.x(), vtxloc.y(), 0);
      os << "v " << vtx.x() << " " << vtx.y() << " " << vtx.z() << "\n";
    }
    // connect them
    os << "f";
    for (size_t i = 1; i <= bounds->vertices().size(); ++i) {
      os << " " << nVtx + i;
    }
    os << "\n";
    nVtx += bounds->vertices().size();
  }
  os.close();
}

using namespace Gaudi;

DECLARE_COMPONENT(GeoSvc)

GeoSvc::GeoSvc(const std::string& name, ISvcLocator* svc)
    : base_class(name, svc)
    //, m_incidentSvc("IncidentSvc", "GeoSvc")
    , m_trackingGeo(nullptr)
    , m_dd4hepGeo(nullptr)
    //, m_geant4geo(0)
    , m_log(msgSvc(), name) {}

GeoSvc::~GeoSvc() {
  if (m_dd4hepGeo) {
    try {
      m_dd4hepGeo->destroyInstance();
      m_dd4hepGeo = 0;
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

  // Genfit
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(
      0., 0., this->centralMagneticField() * 10.0)); // gentfit uses kilo-Gauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  // create a list of all surfaces in the detector:
  dd4hep::rec::SurfaceManager surfMan( *m_dd4hepGeo ) ;
  debug() << " surface manager " << endmsg;
  const auto sM = surfMan.map("tracker") ;
  debug() << " surface map  size: " << sM->size() << endmsg;
  // setup  dd4hep surface map
  //for( dd4hep::rec::SurfaceMap::const_iterator it = sM->begin() ; it != sM->end() ; ++it ){
  for( const auto& [id, s] :   *sM) {
    //dd4hep::rec::Surface* surf = s ;
    m_surfaceMap[ id ] = dynamic_cast<dd4hep::rec::Surface*>(s) ;
    debug() << " surface : " << *s << endmsg;
    m_detPlaneMap[id] = std::shared_ptr<genfit::DetPlane>(
        new genfit::DetPlane({s->origin().x(), s->origin().y(), s->origin().z()}, {s->u().x(), s->u().y(), s->u().z()},
                             {s->v().x(), s->v().y(), s->v().z()}));
  }

  // ACTS
  m_trackingGeo = std::move(Acts::convertDD4hepDetector(m_dd4hepGeo->world(), Acts::Logging::VERBOSE, Acts::equidistant,
                                                        Acts::equidistant, Acts::equidistant));
  if (m_trackingGeo) {
    draw_surfaces(m_trackingGeo, "tracking_geometry.obj");
    debug() << "visiting all the surfaces  " << endmsg;
    m_trackingGeo->visitSurfaces([this](const Acts::Surface* surface) {
      // for now we just require a valid surface
      if (not surface) {
        info() << "no surface??? " << endmsg;
        return;
      }
      auto det_element =
          dynamic_cast<const Acts::DD4hepDetectorElement*>(surface->associatedDetectorElement());

      if (!det_element) {
        error() << "invalid det_element!!! " << endmsg;
        return;
      }
      // more verbose output is lower enum value
      debug() << " det_element->identifier() " << det_element->identifier() << endmsg;
      auto volman  = m_dd4hepGeo->volumeManager();
      auto vol_ctx = volman.lookupContext(det_element->identifier());
      auto vol_id  = vol_ctx->identifier;

      if (msgLevel() <= MSG::DEBUG) {
        auto de  = vol_ctx->element;
        debug() << de.path() << endmsg;
        debug() << de.placementPath() << endmsg;
      }
      this->m_surfaces.insert_or_assign(vol_id, surface);
    });
  }

  m_magneticField = std::make_shared<const Jug::BField::DD4hepBField>(m_dd4hepGeo);
  for (int z : {0, 1000, 2000, 4000}) {
    auto b = m_magneticField->getField({0.0, 0.0, double(z)});
    debug() << "B(z=" << z << " mm) = " << b.transpose()  << " T" << endmsg;
  }

  return StatusCode::SUCCESS;
}

StatusCode GeoSvc::finalize() { return StatusCode::SUCCESS; }

StatusCode GeoSvc::buildDD4HepGeo() {
  // we retrieve the the static instance of the DD4HEP::Geometry
  m_dd4hepGeo = &(dd4hep::Detector::getInstance());
  m_dd4hepGeo->addExtension<IGeoSvc>(this);

  // load geometry
  for (auto& filename : m_xmlFileNames) {
    m_log << MSG::INFO << "loading geometry from file:  '" << filename << "'" << endmsg;
    m_dd4hepGeo->fromCompact(filename);
  }
  m_dd4hepGeo->volumeManager();
  m_dd4hepGeo->apply("DD4hepVolumeManager", 0, 0);
  m_cellid_converter = std::make_shared<const dd4hep::rec::CellIDPositionConverter>(*m_dd4hepGeo);
  return StatusCode::SUCCESS;
}

dd4hep::Detector* GeoSvc::detector() { return (m_dd4hepGeo); }

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
