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
#include "JugBase/IGeoSvc.h"

// ACTS
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"

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

// Create a test context
//#define CHECK_ROTATION_ANGLE(t, a, tolerance)               \
//  {                                                         \
//    Vector3 v = (*t) * Vector3(1, 0, 0);                    \
//    CHECK_CLOSE_ABS(VectorHelpers::phi(v), (a), tolerance); \
//  }

//using SrfVec = std::vector<std::shared_ptr<const Surface>>;
void draw_surfaces(std::shared_ptr<const Acts::TrackingGeometry> trk_geo, const std::string& fname);

class GeoSvc : public extends<Service, IGeoSvc> {
public:
  using VolumeSurfaceMap = std::unordered_map<uint64_t, const Acts::Surface*>;

private:
  /// ACTS Tracking  Geometry
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeo;

  /// Lookup container for hit surfaces that generate smeared hits
  VolumeSurfaceMap m_surfaces;

  /// Pointer to the interface to the DD4hep geometry
  dd4hep::Detector* m_dd4hepgeo;

  std::shared_ptr<const dd4hep::rec::CellIDPositionConverter> m_cellid_converter =
      nullptr; //(*(m_geoSvc->detector()));

  /// XML-files with the detector description
  Gaudi::Property<std::vector<std::string>> m_xmlFileNames{
      this, "detectors", {}, "Detector descriptions XML-files"};

  /// output
  MsgStream m_log;

public:
  GeoSvc(const std::string& name, ISvcLocator* svc);

  virtual ~GeoSvc();

  virtual StatusCode initialize() final;

  virtual StatusCode finalize() final;

  /** Build the dd4hep geometry.
   * This function generates the DD4hep geometry.
   */
  StatusCode buildDD4HepGeo();

  /** Get the top level DetElement.
   *   DD4hep Geometry
   */
  virtual dd4hep::DetElement getDD4HepGeo() override;

  virtual std::shared_ptr<const dd4hep::rec::CellIDPositionConverter>
  cellIDPositionConverter() const
  {
    return m_cellid_converter;
  }

  /** Get the main dd4hep Detector.
   * Returns the pointer to the main dd4hep detector class.
   */
  virtual dd4hep::Detector* detector() override;

  /** Gets the ACTS tracking geometry.
   */
  virtual std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const;

  virtual double centralMagneticField() const
  {
    return m_dd4hepgeo->field().magneticField({0, 0, 0}).z() *
           (Acts::UnitConstants::T / dd4hep::tesla);
  }

  virtual const VolumeSurfaceMap& surfaceMap() const { return m_surfaces; }
};

inline std::shared_ptr<const Acts::TrackingGeometry> GeoSvc::trackingGeometry() const
{
  return m_trackingGeo;
}

#endif // GEOSVC_H
