#include "JugTrack/GeometryContainers.hpp"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"
#include "DD4hep/Volumes.h"
#include "DD4hep/DD4hepUnits.h"


#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"

#include "JugTrack/Index.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Measurement.hpp"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

  /** Source source Linker.
   *
   * The source linker creates "source links" which map the hit to the tracking surface.
   * It also creates "measurements" which take the hit information and creates a corresponding
   * "measurement" which contains the covariance matrix and other geometry related hit information.
   *
   * \ingroup tracking
   */
  class TrackerSourceLinker : public GaudiAlgorithm {
  public:
    DataHandle<eic::TrackerHitCollection>    m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<std::list<IndexSourceLink>>   m_sourceLinkStorage{"sourceLinkStorage", Gaudi::DataHandle::Writer, this};
    DataHandle<IndexSourceLinkContainer>     m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
    DataHandle<MeasurementContainer>         m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

  public:
    TrackerSourceLinker(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("sourceLinkStorage",  m_sourceLinkStorage, "");
      declareProperty("outputSourceLinks",  m_outputSourceLinks, "");
      declareProperty("outputMeasurements", m_outputMeasurements, "");
    }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      constexpr double mm_acts = Acts::UnitConstants::mm;
      constexpr double mm_conv = mm_acts / dd4hep::mm; // = 1/0.1

      // input collection
      const eic::TrackerHitCollection* hits = m_inputHitCollection.get();
      // Create output collections
      auto linkStorage  = m_sourceLinkStorage.createAndPut();
      auto sourceLinks  = m_outputSourceLinks.createAndPut();
      auto measurements = m_outputMeasurements.createAndPut();
      sourceLinks->reserve(hits->size());
      measurements->reserve(hits->size());

      if (msgLevel(MSG::DEBUG)) {
        debug() << (*hits).size() << " hits " << endmsg;
      }
      int ihit = 0;
      for(const auto& ahit : *hits) {

        Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
        cov(0, 0)            = ahit.positionError().xx * mm_acts * mm_acts; // note mm = 1 (Acts)
        cov(1, 1)            = ahit.positionError().yy * mm_acts * mm_acts;
        if (msgLevel(MSG::DEBUG)) {
          debug() << "cov matrix:\n" << cov << endmsg;
        }

        auto vol_ctx   = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
        auto vol_id    = vol_ctx->identifier;
        auto volman    = m_geoSvc->detector()->volumeManager();
        auto alignment = volman.lookupDetElement(vol_id).nominal();
        auto local_position =
            (alignment.worldToLocal({ahit.position().x/mm_conv , ahit.position().y/mm_conv , ahit.position().z/mm_conv }))*mm_conv;

        const auto is      = m_geoSvc->surfaceMap().find(vol_id);
        if (is == m_geoSvc->surfaceMap().end()) {
          error() << " vol_id (" <<  vol_id << ")  not found in m_surfaces." <<endmsg;
          continue;
        }
        const Acts::Surface* surface = is->second;
        // variable surf_center not used anywhere;
        //auto surf_center = surface->center(Acts::GeometryContext());

        // transform global position into local coordinates
        // geometry context contains nothing here
        Acts::Vector2 pos = surface
                                ->globalToLocal(Acts::GeometryContext(),
                                                {ahit.position().x, ahit.position().y, ahit.position().z}, {0, 0, 0})
                                .value();

        Acts::Vector2 loc     = Acts::Vector2::Zero();
        loc[Acts::eBoundLoc0] = pos[0];
        loc[Acts::eBoundLoc1] = pos[1];

        if (msgLevel(MSG::DEBUG)) {
          debug() << " hit position     : " << ahit.position().x << " " << ahit.position().y << " " << ahit.position().z << endmsg;
          debug() << " dd4hep loc pos   : " << local_position.x() << " " << local_position.y() << " " << local_position.z() << endmsg;
          debug() << " surface center   :" << surface->center(Acts::GeometryContext()).transpose() << endmsg;
          debug() << " acts local center:" << pos.transpose() << endmsg;
          debug() << " acts loc pos     : " << loc[Acts::eBoundLoc0] << ", " << loc[Acts::eBoundLoc1] << endmsg;
        }

        // the measurement container is unordered and the index under which the
        // measurement will be stored is known before adding it.
        //
        // variable hitIdx not used anywhere
        //Index hitIdx = measurements->size();
        linkStorage->emplace_back(surface->geometryId(), ihit);
        IndexSourceLink& sourceLink = linkStorage->back();
        auto meas = Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);

        // add to output containers. since the input is already geometry-order,
        // new elements in geometry containers can just be appended at the end.
        sourceLinks->emplace_hint(sourceLinks->end(), sourceLink);
        measurements->emplace_back(std::move(meas));

        ihit++;
      }
      return StatusCode::SUCCESS;
    }

   };
  DECLARE_COMPONENT(TrackerSourceLinker)

} // namespace Jug::reco

