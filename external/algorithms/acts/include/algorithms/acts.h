// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
#pragma once

#include <algorithms/geo.h>
#include <algorithms/logger.h>
#include <algorithms/service.h>

#include "Acts/Definitions/Common.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace algorithms {

class ActsSvc : public LoggedService<ActsSvc> {
public:
  using VolumeSurfaceMap = std::unordered_map<uint64_t, const Acts::Surface*>;

  void init() {
    ;
    ;
  }

private:
  // Getting a handle to the GeoSvc automatically adds the GeoSvc to the Service stack
  // before this ActsSvc. In other words, getting a handle implicitly specifies a [needs]
  // relation.
  const GeoSvc& m_geo = GeoSvc::instance();

  Acts::GeometryContext m_trackingGeoCtx;
  std::unique_ptr<const Acts::TrackingGeometry> m_trackingGeo{nullptr};
  std::unique_ptr<const Acts::IMaterialDecorator> m_materialDeco{nullptr};
  VolumeSurfaceMap m_surfaces;

  ALGORITHMS_DEFINE_LOGGED_SERVICE(ActsSvc)
};

} // namespace algorithms

