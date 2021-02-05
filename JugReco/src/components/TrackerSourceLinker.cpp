#include "JugReco/GeometryContainers.hpp"

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

#include "JugReco/Index.hpp"
#include "JugReco/IndexSourceLink.hpp"
#include "JugReco/Measurement.hpp"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

  class TrackerSourceLinker : public GaudiAlgorithm {
  public:
    DataHandle<eic::TrackerHitCollection>    m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<IndexSourceLinkContainer>     m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
    DataHandle<MeasurementContainer>          m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    /// Lookup container for hit surfaces that generate smeared hits
    std::unordered_map<uint64_t, const Acts::Surface*> m_surfaces;

  public:
    TrackerSourceLinker(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputSourceLinks", m_outputSourceLinks, "");
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
      debug() << "visiting all the surfaces  " << endmsg;
      m_geoSvc->trackingGeometry()->visitSurfaces([this](const Acts::Surface* surface) {
        // for now we just require a valid surface
        if (not surface) {
          return;
        }
        auto det_element = dynamic_cast<const Acts::DD4hepDetectorElement*>(surface->associatedDetectorElement());
        if(!det_element) {
          debug() << "invalid det_element!!! " << endmsg;
          return;
        }
        this->m_surfaces.insert_or_assign(det_element->identifier(), surface);
      });

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const eic::TrackerHitCollection* hits = m_inputHitCollection.get();
      // Create output collections
      auto sourceLinks = m_outputSourceLinks.createAndPut();
      auto measurements = m_outputMeasurements.createAndPut();
      // IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
      // IndexMultimap<Index> hitSimHitsMap;
      sourceLinks->reserve(hits->size());
      measurements->reserve(hits->size());

      debug() << (*hits).size() << " hits " << endmsg;
      int ihit = 0;
      for(const auto& ahit : *hits) {

        Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
        cov(0,0) = ahit.covsym_xx()*Acts::UnitConstants::mm*ahit.covsym_xx()*Acts::UnitConstants::mm;
        cov(1,1) = ahit.covsym_yy()*Acts::UnitConstants::mm*ahit.covsym_yy()*Acts::UnitConstants::mm;

        auto       vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
        auto       vol_id  = vol_ctx->identifier;
        const auto is      = m_surfaces.find(vol_id);
        if (is == m_surfaces.end()) {
          debug() << " vol_id (" <<  vol_id << ")  not found in m_surfaces." <<endmsg;
          continue;
        }
        const Acts::Surface* surface = is->second;
        debug() << " surface center : " << surface->center(Acts::GeometryContext()) << endmsg;
        // transform global position into local coordinates
        Acts::Vector2 pos(0, 0);
        // geometry context contains nothing here
        pos = surface->globalToLocal(Acts::GeometryContext(), {ahit.x(), ahit.y(), ahit.z()}, {0, 0, 0}).value();//, pos);

        Acts::Vector2 loc = Acts::Vector2::Zero();
        loc[Acts::eBoundLoc0]     = pos[0] ;//+ m_cfg.sigmaLoc0 * stdNormal(rng);
        loc[Acts::eBoundLoc1]     = pos[1] ;//+ m_cfg.sigmaLoc1 * stdNormal(rng);
        //debug() << "loc : (" << loc[0] << ", " << loc[1] << ")" << endmsg;

        //local position
        //auto loc = {ahit.x(), ahit.y(), ahit.z()} - vol_ctx->volumePlacement().position()
        //debug() << " hit          : \n" <<  ahit << endmsg;
        //debug() << " cell ID : " << ahit.cellID() << endmsg;
        //debug() << " position : (" <<  ahit.position(0) << ", " <<  ahit.position(1) << ", "<<  ahit.position(2) << ") " << endmsg;
        //debug() << " vol_id       : " <<  vol_id << endmsg;
        //debug() << " placment pos : " << vol_ctx->volumePlacement().position() << endmsg;

        // the measurement container is unordered and the index under which the
        // measurement will be stored is known before adding it.
        Index hitIdx = measurements->size();
        IndexSourceLink sourceLink(vol_id, ihit);
        auto meas = Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);

        // add to output containers. since the input is already geometry-order,
        // new elements in geometry containers can just be appended at the end.
        sourceLinks->emplace_hint(sourceLinks->end(), std::move(sourceLink));
        measurements->emplace_back(std::move(meas));

        ihit++;
      }
      return StatusCode::SUCCESS;
    }

   };
  DECLARE_COMPONENT(TrackerSourceLinker)

} // namespace Jug::reco

