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

#include "JugTrack/Index.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/ProtoTrack.hpp"
#include "JugTrack/GeometryContainers.hpp"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

  /** Single Track source Linker and proto tracks.
   *
   * This algorithm assumes only single track events.
   *
   * \ingroup tracking
   */
  class SingleTrackSourceLinker : public GaudiAlgorithm {
  public:
    DataHandle<eic::TrackerHitCollection>  m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<std::list<IndexSourceLink>> m_sourceLinkStorage{"sourceLinkStorage", Gaudi::DataHandle::Writer, this};
    DataHandle<IndexSourceLinkContainer>   m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
    DataHandle<MeasurementContainer>       m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
    DataHandle<ProtoTrackContainer>        m_outputProtoTracks{"outputProtoTracks", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

  public:
    SingleTrackSourceLinker(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("sourceLinkStorage", m_sourceLinkStorage, "");
      declareProperty("outputSourceLinks", m_outputSourceLinks, "");
      declareProperty("outputMeasurements", m_outputMeasurements, "");
      declareProperty("outputProtoTracks", m_outputProtoTracks, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collection
      const eic::TrackerHitCollection* hits = m_inputHitCollection.get();
      // Create output collections
      auto linkStorage  = m_sourceLinkStorage.createAndPut();
      auto sourceLinks  = m_outputSourceLinks.createAndPut();
      auto measurements = m_outputMeasurements.createAndPut();
      auto protoTracks  = m_outputProtoTracks.createAndPut();
      // IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
      // IndexMultimap<Index> hitSimHitsMap;
      sourceLinks->reserve(hits->size());
      measurements->reserve(hits->size());

      // assume single track --> one ProtoTrack
      ProtoTrack track;
      track.reserve((*hits).size());

      if (msgLevel(MSG::DEBUG)) {
        debug() << (*hits).size() << " hits " << endmsg;
      }
      int ihit = 0;
      for (const auto& ahit : *hits) {

        track.emplace_back(ihit);

        auto       vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
        auto       vol_id  = vol_ctx->identifier;
        const auto is      = m_geoSvc->surfaceMap().find(vol_id);
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

        auto volman    = m_geoSvc->detector()->volumeManager();
        auto detelem   = volman.lookupDetElement(vol_id);
        auto local_pos = detelem.nominal().worldToLocal({ahit.position().x, ahit.position().y, ahit.position().z});
        auto acts_pos  = surface
                            ->globalToLocal(Acts::GeometryContext(),
                                            {ahit.position().x, ahit.position().y, ahit.position().z}, {0, 0, 0})
                            .value();

        if (msgLevel(MSG::DEBUG)) {
          debug() << "===== Debugging hit =====" << endmsg;
          debug() << "DD4hep global pos (" << ahit.position().x << "," << ahit.position().y << "," << ahit.position().z
                  << ")" << endmsg;
          debug() << "DD4hep local  pos (" << local_pos.x() << "," << local_pos.y() << "," << local_pos.z() << ")"
                  << endmsg;
          debug() << "ACTS local position : (" << acts_pos[0] << "," << acts_pos[1] << ")" << endmsg;
          debug() << "ACTS surface center : " << surface->center(Acts::GeometryContext()).transpose() << endmsg;
          debug() << "DD4hep DetElement center : "
                  << detelem.nominal().localToWorld(detelem.placement().position()) / dd4hep::mm << endmsg;
        }
        // construct the vector of measured parameters (2d position in this case)
        Acts::Vector2 pos(acts_pos.x(), acts_pos.y());

        // construct the covariance matrix
        Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
        cov(0, 0) =
            ahit.covMatrix().xx * Acts::UnitConstants::mm * ahit.covMatrix().xx * Acts::UnitConstants::mm;
        cov(1, 1) =
            ahit.covMatrix().yy * Acts::UnitConstants::mm * ahit.covMatrix().yy * Acts::UnitConstants::mm;

        // Above we only consider the two position coordinates the comment below shows how to add time
        // which we will probably want to try later.
        //
        // Acts::SymMatrix3 cov;
        // cov << 0.05, 0., 0., 0., 0.05, 0., 0., 0., 900. * Acts::UnitConstants::ps * Acts::UnitConstants::ps;
        // Acts::Vector3 par(localX, localY, simHit.time());

        linkStorage->emplace_back(surface->geometryId(), ihit);
        IndexSourceLink& sourceLink = linkStorage->back();
        auto            meas =
            Acts::makeMeasurement(sourceLink, pos, cov, Acts::eBoundLoc0, Acts::eBoundLoc1); //, Acts::eBoundTime);

        // TODO: check that local  to global is the same in dd4hep and acts.
        //
        // debug() << "ACTS surface center : " << surface->center(Acts::GeometryContext()) << endmsg;
        // transform global position into local coordinates
        // geometry context contains nothing here
        // Acts::Vector2 loc = Acts::Vector2::Zero();
        // loc[Acts::eBoundLoc0]     = pos[0] ;//+ m_cfg.sigmaLoc0 * stdNormal(rng);
        // loc[Acts::eBoundLoc1]     = pos[1] ;//+ m_cfg.sigmaLoc1 * stdNormal(rng);
        ////debug() << "loc : (" << loc[0] << ", " << loc[1] << ")" << endmsg;

        ////local position
        ////auto loc = {ahit.x(), ahit.y(), ahit.z()} - vol_ctx->volumePlacement().position()
        ////debug() << " hit          : \n" <<  ahit << endmsg;
        ////debug() << " cell ID : " << ahit.cellID() << endmsg;
        ////debug() << " position : (" <<  ahit.position(0) << ", " <<  ahit.position(1) << ", "<<  ahit.position(2) <<
        ///") " << endmsg; /debug() << " vol_id       : " <<  vol_id << endmsg; /debug() << " placment pos : " <<
        /// vol_ctx->volumePlacement().position() << endmsg;

        //// the measurement container is unordered and the index under which the
        //// measurement will be stored is known before adding it.
        // auto meas = Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);

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
  DECLARE_COMPONENT(SingleTrackSourceLinker)

} // namespace Jug::Reco

