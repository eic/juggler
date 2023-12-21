// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#ifndef IACTSGEOSVC_H
#define IACTSGEOSVC_H

#include <GaudiKernel/IService.h>
#include <unordered_map>

namespace dd4hep {
  class Detector;
  class DetElement;
  namespace rec {
    class CellIDPositionConverter;
    class Surface;
  }
} // namespace dd4hep

namespace Acts {
  class TrackingGeometry;
  class Surface;
  class MagneticFieldProvider;
}

/** Geometry service interface.
 *
 * \ingroup base
 * \ingroup geosvc
 */
class GAUDI_API IActsGeoSvc : virtual public IService {
public:
  using VolumeSurfaceMap = std::unordered_map<uint64_t, const Acts::Surface*>;

public:
  /// InterfaceID
  DeclareInterfaceID(IActsGeoSvc, 1, 0);
  virtual std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const = 0;

  virtual std::shared_ptr<const Acts::MagneticFieldProvider>  getFieldProvider() const = 0;
  virtual double centralMagneticField() const = 0;

  virtual const VolumeSurfaceMap& surfaceMap() const = 0;

  virtual ~IActsGeoSvc() {}
};

#endif  // IACTSGEOSVC_H
