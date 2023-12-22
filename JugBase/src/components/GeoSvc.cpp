// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

#include "GeoSvc.h"
#include "GaudiKernel/Service.h"
#include "TGeoManager.h"

#include "DD4hep/Printout.h"

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
  if (buildDD4HepGeo().isFailure()) {
    m_log << MSG::ERROR << "Could not build DD4Hep geometry" << endmsg;
  } else {
    m_log << MSG::INFO << "DD4Hep geometry SUCCESSFULLY built" << endmsg;
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
  return StatusCode::SUCCESS;
}

dd4hep::Detector* GeoSvc::getDetector() { return (m_dd4hepGeo); }

dd4hep::DetElement GeoSvc::getDD4HepGeo() { return (getDetector()->world()); }
