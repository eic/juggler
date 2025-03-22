// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck

#include "ActsExamples/EventData/GeometryContainers.hpp"

// Gaudi
#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>
#include "JugTrack/IActsGeoSvc.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Volumes.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#if Acts_VERSION_MAJOR < 36
#include "Acts/EventData/Measurement.hpp"
#endif
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include "edm4eic/TrackerHitCollection.h"

namespace Jug::Reco {

/** Source source Linker.
 *
 * The source linker creates "source links" which map the hit to the tracking surface.
 * It also creates "measurements" which take the hit information and creates a corresponding
 * "measurement" which contains the covariance matrix and other geometry related hit information.
 *
 * \ingroup tracking
 */
class TrackerSourceLinker : public Gaudi::Algorithm {
private:
  mutable DataHandle<edm4eic::TrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<std::list<ActsExamples::IndexSourceLink>> m_sourceLinkStorage{"sourceLinkStorage", Gaudi::DataHandle::Writer, this};
#if Acts_VERSION_MAJOR < 37 || (Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR == 0)
  mutable DataHandle<ActsExamples::IndexSourceLinkContainer> m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
#endif
  mutable DataHandle<ActsExamples::MeasurementContainer> m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
  /// Pointer to the geometry services
  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IActsGeoSvc> m_actsGeoSvc;
  std::shared_ptr<const dd4hep::rec::CellIDPositionConverter> m_converter;

public:
  TrackerSourceLinker(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("sourceLinkStorage", m_sourceLinkStorage, "");
#if Acts_VERSION_MAJOR < 37 || (Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR == 0)
    declareProperty("outputSourceLinks", m_outputSourceLinks, "");
#endif
    declareProperty("outputMeasurements", m_outputMeasurements, "");
  }

  StatusCode initialize() override {
    if (Gaudi::Algorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc in the right place in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    m_actsGeoSvc = service("ActsGeoSvc");
    if (!m_actsGeoSvc) {
      error() << "Unable to locate ACTS Geometry Service. "
              << "Make sure you have ActsGeoSvc in the right place in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    m_converter = std::make_shared<const dd4hep::rec::CellIDPositionConverter>(*(m_geoSvc->getDetector()));

    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    constexpr double mm_acts = Acts::UnitConstants::mm;
    constexpr double mm_conv = mm_acts / dd4hep::mm; // = 1/0.1

    // input collection
    const edm4eic::TrackerHitCollection* hits = m_inputHitCollection.get();
    // Create output collections
    auto* linkStorage  = m_sourceLinkStorage.createAndPut();
#if Acts_VERSION_MAJOR < 37 || (Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR == 0)
    auto* sourceLinks  = m_outputSourceLinks.createAndPut();
    sourceLinks->reserve(hits->size());
#endif
    auto* measurements = m_outputMeasurements.createAndPut();
    measurements->reserve(hits->size());

    if (msgLevel(MSG::DEBUG)) {
      debug() << (*hits).size() << " hits " << endmsg;
    }
    int ihit = 0;
    for (const auto& ahit : *hits) {

      Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Zero();
      cov(0, 0)            = ahit.getPositionError().xx * mm_acts * mm_acts; // note mm = 1 (Acts)
      cov(1, 1)            = ahit.getPositionError().yy * mm_acts * mm_acts;
      if (msgLevel(MSG::DEBUG)) {
        debug() << "cov matrix:\n" << cov << endmsg;
      }

      const auto* vol_ctx = m_converter->findContext(ahit.getCellID());
      auto vol_id = vol_ctx->identifier;

      const auto is = m_actsGeoSvc->surfaceMap().find(vol_id);
      if (is == m_actsGeoSvc->surfaceMap().end()) {
        error() << " vol_id (" << vol_id << ")  not found in m_surfaces." << endmsg;
        continue;
      }
      const Acts::Surface* surface = is->second;
      // variable surf_center not used anywhere;
      // auto surf_center = surface->center(Acts::GeometryContext());

      // transform global position into local coordinates
      // geometry context contains nothing here
      Acts::Vector2 pos =
          surface
              ->globalToLocal(Acts::GeometryContext(),
                              {ahit.getPosition().x, ahit.getPosition().y, ahit.getPosition().z}, {0, 0, 0})
              .value();

      Acts::Vector2 loc     = Acts::Vector2::Zero();
      loc[Acts::eBoundLoc0] = pos[0];
      loc[Acts::eBoundLoc1] = pos[1];

      if (msgLevel(MSG::DEBUG)) {
        auto volman         = m_geoSvc->getDetector()->volumeManager();
        auto alignment      = volman.lookupDetElement(vol_id).nominal();
        auto local_position = (alignment.worldToLocal({ahit.getPosition().x / mm_conv, ahit.getPosition().y / mm_conv,
                                                       ahit.getPosition().z / mm_conv})) *
                              mm_conv;
        debug() << " hit position     : " << ahit.getPosition().x << " " << ahit.getPosition().y << " "
                << ahit.getPosition().z << endmsg;
        debug() << " dd4hep loc pos   : " << local_position.x() << " " << local_position.y() << " "
                << local_position.z() << endmsg;
        debug() << " surface center   :" << surface->center(Acts::GeometryContext()).transpose() << endmsg;
        debug() << " acts local center:" << pos.transpose() << endmsg;
        debug() << " acts loc pos     : " << loc[Acts::eBoundLoc0] << ", " << loc[Acts::eBoundLoc1] << endmsg;
      }

      // the measurement container is unordered and the index under which the
      // measurement will be stored is known before adding it.
      //
      // variable hitIdx not used anywhere
      // Index hitIdx = measurements->size();

      auto geoId = surface->geometryId();

      linkStorage->emplace_back(geoId, ihit);
#if Acts_VERSION_MAJOR < 37 || (Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR == 0)
      ActsExamples::IndexSourceLink& sourceLink = linkStorage->back();
      sourceLinks->emplace_hint(sourceLinks->end(), sourceLink);
#endif

#if Acts_VERSION_MAJOR > 37 || (Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR >= 1)
      std::array<Acts::BoundIndices, 2> indices = {Acts::eBoundLoc0, Acts::eBoundLoc1};
      Acts::visit_measurement(
        indices.size(), [&](auto dim) -> ActsExamples::VariableBoundMeasurementProxy {
          if constexpr (dim == indices.size()) {
            return ActsExamples::VariableBoundMeasurementProxy{
              measurements->emplaceMeasurement<dim>(geoId, indices, loc, cov)
            };
          } else {
            throw std::runtime_error("Dimension not supported in measurement creation");
          }
        }
      );
#elif Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR == 0
      std::array<Acts::BoundIndices, 2> indices = {Acts::eBoundLoc0, Acts::eBoundLoc1};
      Acts::visit_measurement(
        indices.size(), [&](auto dim) -> ActsExamples::VariableBoundMeasurementProxy {
          if constexpr (dim == indices.size()) {
            return ActsExamples::VariableBoundMeasurementProxy{
              measurements->emplaceMeasurement<dim>(Acts::SourceLink{sourceLink}, indices, loc, cov)
            };
          } else {
            throw std::runtime_error("Dimension not supported in measurement creation");
          }
        }
      );
#elif Acts_VERSION_MAJOR == 36 && Acts_VERSION_MINOR >= 1
      auto measurement = ActsExamples::makeVariableSizeMeasurement(
        Acts::SourceLink{sourceLink}, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
      measurements->emplace_back(std::move(measurement));
#elif Acts_VERSION_MAJOR == 36 && Acts_VERSION_MINOR == 0
      auto measurement = ActsExamples::makeFixedSizeMeasurement(
        Acts::SourceLink{sourceLink}, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
      measurements->emplace_back(std::move(measurement));
#else
      auto measurement = Acts::makeMeasurement(
        Acts::SourceLink{sourceLink}, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
      measurements->emplace_back(std::move(measurement));
#endif

      ihit++;
    }
    return StatusCode::SUCCESS;
  }
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(TrackerSourceLinker)

} // namespace Jug::Reco
