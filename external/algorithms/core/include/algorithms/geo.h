// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// DD4hep Geometry service exposing a detector(), world(), and cellIDPositionConverter()
// Meant to be set by the calling framework, but can also load DD4hep itself
// when given a compact file as Property.
//
#pragma once

#include <gsl/gsl>
#include <memory>

#include "DDRec/CellIDPositionConverter.h"
#include <DD4hep/Detector.h>

#include <algorithms/logger.h>
#include <algorithms/service.h>

namespace algorithms {

class GeoSvc : public LoggedService<GeoSvc> {
public:
  // Initialize the geometry service, to be called after with a global detector
  // pointer (will not re-initialize), or after at least the detectors property was specified
  // (will load XML and then fully initialize)
  void init(dd4hep::Detector* = nullptr);

  // TODO check const-ness
  gsl::not_null<const dd4hep::Detector*> detector() const { return m_detector; }
  dd4hep::DetElement world() const { return detector()->world(); }
  gsl::not_null<const dd4hep::rec::CellIDPositionConverter*> cellIDPositionConverter() const {
    return m_converter.get();
  }

private:
  dd4hep::Detector* m_detector = nullptr;
  std::unique_ptr<const dd4hep::rec::CellIDPositionConverter> m_converter;

  // Configuration variables. These only need to be specified if we are actually running
  // in "standalone" mode (hence they're optional)
  Property<std::vector<std::string>> m_xml_list{this, "detectors", {}};

  ALGORITHMS_DEFINE_LOGGED_SERVICE(GeoSvc)
};

} // namespace algorithms
