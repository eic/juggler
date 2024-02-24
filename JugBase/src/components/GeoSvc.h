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
#include <k4Interface/IGeoSvc.h>

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



class GeoSvc : public extends<Service, IGeoSvc> {
private:

  /** DD4hep detector interface class.
   * This is the main dd4hep detector handle.
   * <a href="https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html">See DD4hep Detector documentation</a>
   */
  dd4hep::Detector* m_dd4hepGeo = nullptr;

  /// XML-files with the detector description
  Gaudi::Property<std::vector<std::string>> m_xmlFileNames{
      this, "detectors", {}, "Detector descriptions XML-files"};

  Gaudi::Property<std::string> _{
      this, "materials", "", "(unused)"};

  /// Gaudi logging output
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

  /** Get the main dd4hep Detector.
   * Returns the pointer to the main dd4hep detector class.
   * <a href="https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html">See DD4hep Detector documentation</a>
   */
  virtual dd4hep::Detector* getDetector() override;
  [[deprecated("Use getDetector() instead")]]
  virtual dd4hep::Detector* lcdd() {
    return getDetector();
  };
  virtual G4VUserDetectorConstruction* getGeant4Geo() override {
    assert(false && "getGeant4Geo() noy implemented");
    return nullptr; 
  };
  virtual std::string constantAsString(std::string const& name) override {
    assert(false && "constantsAsString() noy implemented");
    (void)name;
    return "";
  };

};

#endif // GEOSVC_H
