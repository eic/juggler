#include "GeoSvc.h"
#include "GaudiKernel/Service.h"
//#include "GeoConstruction.h"
#include "TGeoManager.h"

#include "DD4hep/Printout.h"

#include "JugBase/ACTSLogger.h"
#include "JugBase/Acts/MaterialWiper.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
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

static const std::map<int, Acts::Logging::Level> s_msgMap = {
    {MSG::DEBUG, Acts::Logging::DEBUG},
    {MSG::VERBOSE, Acts::Logging::VERBOSE},
    {MSG::INFO, Acts::Logging::INFO},
    {MSG::WARNING, Acts::Logging::WARNING},
    {MSG::FATAL, Acts::Logging::FATAL},
    {MSG::ERROR, Acts::Logging::ERROR},
};

void draw_surfaces(std::shared_ptr<const Acts::TrackingGeometry> trk_geo, const Acts::GeometryContext geo_ctx, const std::string& fname)
{
  using namespace Acts;
  std::vector<const Surface*> surfaces;

  trk_geo->visitSurfaces([&](const Acts::Surface* surface) {
    // for now we just require a valid surface
    if (surface == nullptr) {
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
    const auto* srf    = dynamic_cast<const PlaneSurface*>(srfx);
    const auto* bounds = dynamic_cast<const PlanarBounds*>(&srf->bounds());
    for (const auto& vtxloc : bounds->vertices()) {
      Vector3 vtx = srf->transform(geo_ctx) * Vector3(vtxloc.x(), vtxloc.y(), 0);
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

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(GeoSvc)

GeoSvc::GeoSvc(const std::string& name, ISvcLocator* svc)
    : base_class(name, svc)
    , m_log(msgSvc(), name) {}

GeoSvc::~GeoSvc() {
  if (m_dd4hepGeo != nullptr) {
    try {
      m_dd4hepGeo->destroyInstance();
      m_dd4hepGeo = nullptr;
    } catch (...) {
    }
  }
}

StatusCode GeoSvc::initialize() {
  StatusCode sc = Service::initialize();
  if (!sc.isSuccess()) {
    return sc;
  }
  // Turn off TGeo printouts if appropriate for the msg level
  if (msgLevel() >= MSG::INFO) {
    TGeoManager::SetVerboseLevel(0);
  }
  uint printoutLevel = msgLevel();
  dd4hep::setPrintLevel(dd4hep::PrintLevel(printoutLevel));
  // m_incidentSvc->addListener(this, "GeometryFailure");
  if (buildDD4HepGeo().isFailure()) {
    m_log << MSG::ERROR << "Could not build DD4Hep geometry" << endmsg;
  } else {
    m_log << MSG::INFO << "DD4Hep geometry SUCCESSFULLY built" << endmsg;
  }

  // Genfit
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(
      0., 0., this->centralMagneticField() * 10.0)); // gentfit uses kilo-Gauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  // create a list of all surfaces in the detector:
  dd4hep::rec::SurfaceManager surfMan( *m_dd4hepGeo ) ;
  debug() << " surface manager " << endmsg;
  const auto* const sM = surfMan.map("tracker") ;
  if (sM != nullptr) {
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
  }

  // Set ACTS logging level
  auto im = s_msgMap.find(msgLevel());
  if (im != s_msgMap.end()) {
    m_actsLoggingLevel = im->second;
  }

  // Load ACTS materials maps
  if (m_jsonFileName.size() > 0) {
    m_log << MSG::INFO << "loading materials map from file:  '" << m_jsonFileName << "'" << endmsg;
    // Set up the converter first
    Acts::MaterialMapJsonConverter::Config jsonGeoConvConfig;
    // Set up the json-based decorator
    m_materialDeco = std::make_shared<const Acts::JsonMaterialDecorator>(
      jsonGeoConvConfig, m_jsonFileName, m_actsLoggingLevel);
  } else {
    m_log << MSG::WARNING << "no ACTS materials map has been loaded" << endmsg;
    m_materialDeco = std::make_shared<const Acts::MaterialWiper>();
  }

  // Convert DD4hep geometry to ACTS
  Acts::BinningType bTypePhi = Acts::equidistant;
  Acts::BinningType bTypeR = Acts::equidistant;
  Acts::BinningType bTypeZ = Acts::equidistant;
  double layerEnvelopeR = Acts::UnitConstants::mm;
  double layerEnvelopeZ = Acts::UnitConstants::mm;
  double defaultLayerThickness = Acts::UnitConstants::fm;
  using Acts::sortDetElementsByID;
  m_trackingGeo = Acts::convertDD4hepDetector(
      m_dd4hepGeo->world(),
      m_actsLoggingLevel,
      bTypePhi,
      bTypeR,
      bTypeZ,
      layerEnvelopeR,
      layerEnvelopeZ,
      defaultLayerThickness,
      sortDetElementsByID,
      m_trackingGeoCtx,
      m_materialDeco);
  // Visit surfaces
  if (m_trackingGeo) {
    draw_surfaces(m_trackingGeo, m_trackingGeoCtx, "tracking_geometry.obj");
    debug() << "visiting all the surfaces  " << endmsg;
    m_trackingGeo->visitSurfaces([this](const Acts::Surface* surface) {
      // for now we just require a valid surface
      if (surface == nullptr) {
        info() << "no surface??? " << endmsg;
        return;
      }
      const auto* det_element =
        dynamic_cast<const Acts::DD4hepDetectorElement*>(surface->associatedDetectorElement());

      if (det_element == nullptr) {
        error() << "invalid det_element!!! " << endmsg;
        return;
      }
      // more verbose output is lower enum value
      debug() << " det_element->identifier() " << det_element->identifier() << endmsg;
      auto volman  = m_dd4hepGeo->volumeManager();
      auto* vol_ctx = volman.lookupContext(det_element->identifier());
      auto vol_id  = vol_ctx->identifier;

      if (msgLevel() <= MSG::DEBUG) {
        auto de  = vol_ctx->element;
        debug() << de.path() << endmsg;
        debug() << de.placementPath() << endmsg;
      }
      this->m_surfaces.insert_or_assign(vol_id, surface);
    });
  }

  // Load ACTS magnetic field
  m_magneticField = std::make_shared<const Jug::BField::DD4hepBField>(m_dd4hepGeo);
  Acts::MagneticFieldContext m_fieldctx{Jug::BField::BFieldVariant(m_magneticField)};
  auto bCache = m_magneticField->makeCache(m_fieldctx);
  for (int z : {0, 1000, 2000, 4000}) {
    auto b = m_magneticField->getField({0.0, 0.0, double(z)}, bCache).value();
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
  m_dd4hepGeo->apply("DD4hepVolumeManager", 0, nullptr);
  m_cellid_converter = std::make_shared<const dd4hep::rec::CellIDPositionConverter>(*m_dd4hepGeo);
  return StatusCode::SUCCESS;
}

dd4hep::Detector* GeoSvc::detector() { return (m_dd4hepGeo); }

dd4hep::DetElement GeoSvc::getDD4HepGeo() { return (detector()->world()); }
