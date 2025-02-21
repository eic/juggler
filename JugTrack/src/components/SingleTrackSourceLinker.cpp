// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck

// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
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

#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"

#include "edm4eic/TrackerHitCollection.h"

namespace Jug::Reco {

/** Single Track source Linker and proto tracks.
 *
 * This algorithm assumes only single track events.
 *
 * \ingroup tracking
 */
class SingleTrackSourceLinker : public GaudiAlgorithm {
private:
  DataHandle<edm4eic::TrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<std::list<ActsExamples::IndexSourceLink>> m_sourceLinkStorage{"sourceLinkStorage", Gaudi::DataHandle::Writer, this};
#if Acts_VERSION_MAJOR < 37 || (Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR == 0)
  DataHandle<ActsExamples::IndexSourceLinkContainer> m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
#endif
  DataHandle<ActsExamples::MeasurementContainer> m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
  DataHandle<ActsExamples::ProtoTrackContainer> m_outputProtoTracks{"outputProtoTracks", Gaudi::DataHandle::Writer, this};
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IActsGeoSvc> m_actsGeoSvc;
  std::shared_ptr<const dd4hep::rec::CellIDPositionConverter> m_converter;

public:
  SingleTrackSourceLinker(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("sourceLinkStorage", m_sourceLinkStorage, "");
#if Acts_VERSION_MAJOR < 37 || (Acts_VERSION_MAJOR == 37 && Acts_VERSION_MINOR == 0)
    declareProperty("outputSourceLinks", m_outputSourceLinks, "");
#endif
    declareProperty("outputMeasurements", m_outputMeasurements, "");
    declareProperty("outputProtoTracks", m_outputProtoTracks, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    m_actsGeoSvc = service("ActsGeoSvc");
    if (!m_actsGeoSvc) {
      error() << "Unable to locate ACTS Geometry Service. "
              << "Make sure you have ActsGeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    m_converter = std::make_shared<const dd4hep::rec::CellIDPositionConverter>(*(m_geoSvc->getDetector()));
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
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
    auto* protoTracks  = m_outputProtoTracks.createAndPut();
    // IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
    // IndexMultimap<Index> hitSimHitsMap;

    // assume single track --> one ProtoTrack
    ActsExamples::ProtoTrack track;
    track.reserve((*hits).size());

    if (msgLevel(MSG::DEBUG)) {
      debug() << (*hits).size() << " hits " << endmsg;
    }
    int ihit = 0;
    for (const auto& ahit : *hits) {

      track.emplace_back(ihit);

      const auto* vol_ctx = m_converter->findContext(ahit.getCellID());
      auto vol_id   = vol_ctx->identifier;
      const auto is = m_actsGeoSvc->surfaceMap().find(vol_id);
      if (is == m_actsGeoSvc->surfaceMap().end()) {
        if (msgLevel(MSG::DEBUG)) {
          debug() << " vol_id (" << vol_id << ")  not found in m_surfaces!!!!" << endmsg;
        }
        continue;
      }
      const Acts::Surface* surface = is->second;

      // NOTE
      // Here is where the important hit and tracking geometry is connected.
      // A "Measurement" is constructed to for each hit which makes the connection to
      // the tracking surface and covariance matrix

      auto acts_pos = surface
                          ->globalToLocal(Acts::GeometryContext(),
                                          {ahit.getPosition().x, ahit.getPosition().y, ahit.getPosition().z}, {0, 0, 0})
                          .value();

      if (msgLevel(MSG::DEBUG)) {
        auto volman  = m_geoSvc->getDetector()->volumeManager();
        auto detelem = volman.lookupDetElement(vol_id);
        auto local_pos =
            detelem.nominal().worldToLocal({ahit.getPosition().x, ahit.getPosition().y, ahit.getPosition().z});
        debug() << "===== Debugging hit =====" << endmsg;
        debug() << "DD4hep global pos (" << ahit.getPosition().x << "," << ahit.getPosition().y << ","
                << ahit.getPosition().z << ")" << endmsg;
        debug() << "DD4hep local  pos (" << local_pos.x() << "," << local_pos.y() << "," << local_pos.z() << ")"
                << endmsg;
        debug() << "ACTS local position : (" << acts_pos[0] << "," << acts_pos[1] << ")" << endmsg;
        debug() << "ACTS surface center : " << surface->center(Acts::GeometryContext()).transpose() << endmsg;
        debug() << "DD4hep DetElement center : "
                << detelem.nominal().localToWorld(detelem.placement().position()) / dd4hep::mm << endmsg;
      }
      // construct the vector of measured parameters (2d getPosition in this case)
      Acts::Vector2 pos(acts_pos.x(), acts_pos.y());

      // construct the covariance matrix
      Acts::SquareMatrix2 cov = Acts::SquareMatrix2::Zero();
      cov(0, 0)            = ahit.getPositionError().xx * Acts::UnitConstants::mm * Acts::UnitConstants::mm;
      cov(1, 1)            = ahit.getPositionError().yy * Acts::UnitConstants::mm * Acts::UnitConstants::mm;

      // Above we only consider the two position coordinates the comment below shows how to add time
      // which we will probably want to try later.
      //
      // Acts::SquareMatrix3 cov;
      // cov << 0.05, 0., 0., 0., 0.05, 0., 0., 0., 900. * Acts::UnitConstants::ps * Acts::UnitConstants::ps;
      // Acts::Vector3 par(localX, localY, simHit.time());

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
              measurements->emplaceMeasurement<dim>(geoId, indices, pos, cov)
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
              measurements->emplaceMeasurement<dim>(Acts::SourceLink{sourceLink}, indices, pos, cov)
            };
          } else {
            throw std::runtime_error("Dimension not supported in measurement creation");
          }
        }
      );
#elif Acts_VERSION_MAJOR == 36 && Acts_VERSION_MINOR >= 1
      auto measurement = ActsExamples::makeVariableSizeMeasurement(
          Acts::SourceLink{sourceLink}, pos, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
      measurements->emplace_back(std::move(measurement));
#elif Acts_VERSION_MAJOR == 36 && Acts_VERSION_MINOR == 0
      auto measurement = ActsExamples::makeFixedSizeMeasurement(
          Acts::SourceLink{sourceLink}, pos, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
      measurements->emplace_back(std::move(measurement));
#else
      auto measurement =
          Acts::makeMeasurement(Acts::SourceLink{sourceLink}, pos, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
      measurements->emplace_back(std::move(measurement));
#endif

      ihit++;
    }
    // add proto track to the output collection
    protoTracks->emplace_back(std::move(track));
    return StatusCode::SUCCESS;
  }
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(SingleTrackSourceLinker)

} // namespace Jug::Reco
