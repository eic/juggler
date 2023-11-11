// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck

// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Volumes.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
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
  DataHandle<ActsExamples::IndexSourceLinkContainer> m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
  DataHandle<ActsExamples::MeasurementContainer> m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
  DataHandle<ActsExamples::ProtoTrackContainer> m_outputProtoTracks{"outputProtoTracks", Gaudi::DataHandle::Writer, this};
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

public:
  SingleTrackSourceLinker(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("sourceLinkStorage", m_sourceLinkStorage, "");
    declareProperty("outputSourceLinks", m_outputSourceLinks, "");
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
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collection
    const edm4eic::TrackerHitCollection* hits = m_inputHitCollection.get();
    // Create output collections
    auto* linkStorage  = m_sourceLinkStorage.createAndPut();
    auto* sourceLinks  = m_outputSourceLinks.createAndPut();
    auto* measurements = m_outputMeasurements.createAndPut();
    auto* protoTracks  = m_outputProtoTracks.createAndPut();
    // IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
    // IndexMultimap<Index> hitSimHitsMap;
    sourceLinks->reserve(hits->size());
    measurements->reserve(hits->size());

    // assume single track --> one ProtoTrack
    ActsExamples::ProtoTrack track;
    track.reserve((*hits).size());

    if (msgLevel(MSG::DEBUG)) {
      debug() << (*hits).size() << " hits " << endmsg;
    }
    int ihit = 0;
    for (const auto& ahit : *hits) {

      track.emplace_back(ihit);

      const auto* vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.getCellID());
      auto vol_id   = vol_ctx->identifier;
      const auto is = m_geoSvc->surfaceMap().find(vol_id);
      if (is == m_geoSvc->surfaceMap().end()) {
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
        auto volman  = m_geoSvc->detector()->volumeManager();
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

      linkStorage->emplace_back(surface->geometryId(), ihit);
      ActsExamples::IndexSourceLink& sourceLink = linkStorage->back();
      auto meas =
          Acts::makeMeasurement(Acts::SourceLink{sourceLink}, pos, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);

      // add to output containers. since the input is already geometry-order,
      // new elements in geometry containers can just be appended at the end.
      sourceLinks->emplace_hint(sourceLinks->end(), sourceLink);
      measurements->emplace_back(std::move(meas));

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
