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

// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"

// DD4Hep
#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

class GeoSvc : public extends<Service, IGeoSvc> {

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

  virtual std::shared_ptr<const dd4hep::rec::CellIDPositionConverter> cellIDPositionConverter() const {
    return m_cellid_converter;
  }

  /** Get the main dd4hep Detector.
   * Returns the pointer to the main dd4hep detector class.
   */
  virtual dd4hep::Detector* detector() override ;

  /** Gets the ACTS tracking geometry.
   */
  virtual std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const;

  virtual double centralMagneticField() const  {
    return m_dd4hepgeo->field().magneticField({0,0,0}).z();                                                                                  
  }
 
;


private:

  /// ACTS Tracking  Geometry
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeo;

  /// Pointer to the interface to the DD4hep geometry
  dd4hep::Detector* m_dd4hepgeo;


  std::shared_ptr<const dd4hep::rec::CellIDPositionConverter> m_cellid_converter = nullptr;//(*(m_geoSvc->detector()));
      

  /// XML-files with the detector description
  Gaudi::Property<std::vector<std::string>> m_xmlFileNames{this, "detectors", {}, "Detector descriptions XML-files"};

  /// output
  MsgStream m_log;
};

inline std::shared_ptr<const Acts::TrackingGeometry> GeoSvc::trackingGeometry() const { return m_trackingGeo; }

#endif  // GEOSVC_H
