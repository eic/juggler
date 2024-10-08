// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

#include "ActsGeoSvc.h"
#include "GaudiKernel/Service.h"
#include "TGeoManager.h"

#include "DD4hep/Printout.h"

#include "JugTrack/ACTSLogger.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

static const std::map<int, Acts::Logging::Level> s_msgMap = {
    {MSG::DEBUG, Acts::Logging::DEBUG},
    {MSG::VERBOSE, Acts::Logging::VERBOSE},
    {MSG::INFO, Acts::Logging::INFO},
    {MSG::WARNING, Acts::Logging::WARNING},
    {MSG::ERROR, Acts::Logging::ERROR},
    {MSG::FATAL, Acts::Logging::FATAL},
    {MSG::ALWAYS, Acts::Logging::MAX},
};

using namespace Gaudi;

namespace Jug::Reco {

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ActsGeoSvc)

ActsGeoSvc::ActsGeoSvc(const std::string& name, ISvcLocator* svc)
    : base_class(name, svc)
    , m_log(msgSvc(), name) {}

ActsGeoSvc::~ActsGeoSvc() {
  if (m_dd4hepGeo != nullptr) {
    try {
      m_dd4hepGeo->destroyInstance();
      m_dd4hepGeo = nullptr;
    } catch (...) {
    }
  }
}

StatusCode ActsGeoSvc::initialize() {
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
  }

  // Convert DD4hep geometry to ACTS
  auto logger = Acts::getDefaultLogger("CONV", m_actsLoggingLevel);
  Acts::BinningType bTypePhi = Acts::equidistant;
  Acts::BinningType bTypeR = Acts::equidistant;
  Acts::BinningType bTypeZ = Acts::equidistant;
  double layerEnvelopeR = Acts::UnitConstants::mm;
  double layerEnvelopeZ = Acts::UnitConstants::mm;
  double defaultLayerThickness = Acts::UnitConstants::fm;
  using Acts::sortDetElementsByID;
  m_trackingGeo = Acts::convertDD4hepDetector(
      m_dd4hepGeo->world(),
      *logger,
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

StatusCode ActsGeoSvc::finalize() { return StatusCode::SUCCESS; }

StatusCode ActsGeoSvc::buildDD4HepGeo() {
  // we retrieve the the static instance of the DD4HEP::Geometry
  m_dd4hepGeo = &(dd4hep::Detector::getInstance());
  m_dd4hepGeo->addExtension<IActsGeoSvc>(this);

  // load geometry
  for (auto& filename : m_xmlFileNames) {
    m_log << MSG::INFO << "loading geometry from file:  '" << filename << "'" << endmsg;
    m_dd4hepGeo->fromCompact(filename);
  }
  m_dd4hepGeo->volumeManager();
  m_dd4hepGeo->apply("DD4hepVolumeManager", 0, nullptr);
  return StatusCode::SUCCESS;
}

}
