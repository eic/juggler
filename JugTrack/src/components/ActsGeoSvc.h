// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

//
//  GeoSvc.h
//
//
//  Created by Julia Hrdinka on 30/03/15.
//
//

#ifndef GEOSVC_H
#define GEOSVC_H

// Interface
#include "JugTrack/IActsGeoSvc.h"

// ACTS
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include <Acts/Material/IMaterialDecorator.hpp>

#if Acts_VERSION_MAJOR >= 44
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#else
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#endif

// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"

// DD4Hep
#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"
#include "DD4hep/DD4hepUnits.h"

#include "JugTrack/DD4hepBField.h"


namespace Jug::Reco {

class ActsGeoSvc : public extends<Service, IActsGeoSvc> {
public:
  using VolumeSurfaceMap = std::unordered_map<uint64_t, const Acts::Surface*>;

private:

  /** DD4hep detector interface class.
   * This is the main dd4hep detector handle.
   * <a href="https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html">See DD4hep Detector documentation</a>
   */
  dd4hep::Detector* m_dd4hepGeo = nullptr;

  /// ACTS Logging Level
  Acts::Logging::Level m_actsLoggingLevel = Acts::Logging::INFO;

  /// ACTS Tracking Geometry Context
  Acts::GeometryContext m_trackingGeoCtx;

  /// ACTS Tracking Geometry
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeo{nullptr};

  /// ACTS Material Decorator
  std::shared_ptr<const Acts::IMaterialDecorator> m_materialDeco{nullptr};

  /// ACTS surface lookup container for hit surfaces that generate smeared hits
  VolumeSurfaceMap m_surfaces;

  /// Acts magnetic field
  std::shared_ptr<const Jug::BField::DD4hepBField> m_magneticField = nullptr;

  /// XML-files with the detector description
  Gaudi::Property<std::vector<std::string>> m_xmlFileNames{
      this, "detectors", {}, "Detector descriptions XML-files"};

  /// JSON-file with the material map
  Gaudi::Property<std::string> m_jsonFileName{
      this, "materials", "", "Material map JSON-file"};

  /// Gaudi logging output
  MsgStream m_log;

public:
  ActsGeoSvc(const std::string& name, ISvcLocator* svc);

  virtual ~ActsGeoSvc();

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  /** Build the dd4hep geometry.
   * This function generates the DD4hep geometry.
   */
  StatusCode buildDD4HepGeo();

  /** Gets the ACTS tracking geometry.
   */
  virtual std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const;

  virtual std::shared_ptr<const Acts::MagneticFieldProvider> getFieldProvider() const { return m_magneticField; }

  virtual double centralMagneticField() const
  {
    return m_dd4hepGeo->field().magneticField({0, 0, 0}).z() * (Acts::UnitConstants::T / dd4hep::tesla);
  }

  virtual const VolumeSurfaceMap& surfaceMap() const { return m_surfaces; }
};

inline std::shared_ptr<const Acts::TrackingGeometry> ActsGeoSvc::trackingGeometry() const
{
  return m_trackingGeo;
}

}

#endif // GEOSVC_H
