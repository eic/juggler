// This file is part of the Acts project.
#include "JugTrack/GeometryContainers.hpp"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"
#include "DD4hep/Volumes.h"
#include "DD4hep/DD4hepUnits.h"


#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"

#include "JugTrack/SourceLinks.h"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

  /** Source Linking .
   *
   * \ingroup track
   */
  class TrackingHitsSourceLinker : public GaudiAlgorithm {
  public:
    using TrkHits = eic::TrackerHitCollection;

    Gaudi::Property<std::vector<std::string>> m_inputs{this, "inputTrackerCollections", {}, "List of intput tracker hit collections"};
    std::vector<DataHandle<TrkHits>*>         m_inputHandles;

    // Gaudi::Property<std::vector<std::string>> m_inputCollections{this, "inputCollections", {}, "Input Tracker hit
    // collections"}; DataHandle<eic::TrackerHitCollection>    m_inputHitCollection{"inputHitCollection",
    // Gaudi::DataHandle::Reader, this};
    DataHandle<SourceLinkContainer> m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};

    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    /// Lookup container for hit surfaces that generate smeared hits
    std::unordered_map<uint64_t, const Acts::Surface*> m_surfaces;

  public:
    TrackingHitsSourceLinker(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      // declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputSourceLinks", m_outputSourceLinks, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;

      for (const auto& col : m_inputs) {
        m_inputHandles.push_back(new DataHandle<TrkHits>(col, Gaudi::DataHandle::Reader, this));
        // declareProperty("inputHitCollection", m_inputHitCollection, "");
        declareProperty("dummy_in_" + col, *(m_inputHandles.back()),"");
      }

      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
        return StatusCode::FAILURE;
      }
      debug() << "visiting all the surfaces  " << endmsg;
      m_geoSvc->trackingGeometry()->visitSurfaces([this](const Acts::Surface* surface) {
        // for now we just require a valid surface
        if (not surface) {
          return;
        }
        auto det_element = dynamic_cast<const Acts::DD4hepDetectorElement*>(surface->associatedDetectorElement());
        if (!det_element) {
          debug() << "invalid det_element!!! " << endmsg;
          return;
        }
        this->m_surfaces.insert_or_assign(det_element->identifier(), surface);
      });

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // Create output collections
      auto source_links = m_outputSourceLinks.createAndPut();

      for (const auto& handle : m_inputHandles) {
        // input collection
        const eic::TrackerHitCollection* hits = handle->get();
        debug() << (*hits).size() << " hits " << endmsg;

        for (const auto& ahit : *hits) {

          // setup local covariance
          Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
          cov(Acts::eBoundLoc0, Acts::eBoundLoc0) =
              ahit.xx() * Acts::UnitConstants::mm * ahit.xx() * Acts::UnitConstants::mm;
          cov(Acts::eBoundLoc1, Acts::eBoundLoc1) =
              ahit.yy() * Acts::UnitConstants::mm * ahit.yy() * Acts::UnitConstants::mm;

          auto vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
          auto vol_id  = vol_ctx->identifier;
          debug() << " vol_id       : " << vol_id << endmsg;
          debug() << " placment pos : " << vol_ctx->volumePlacement().position() << endmsg;

          const auto is = m_surfaces.find(vol_id);
          if (is == m_surfaces.end()) {
            debug() << " vol_id (" << vol_id << ")  not found in m_surfaces." << endmsg;
            continue;
          }
          const Acts::Surface* surface = is->second;

          debug() << " surface center : " << surface->center(Acts::GeometryContext()) << endmsg;
          // transform global position into local coordinates
          Acts::Vector2D pos(0, 0);
          // geometry context contains nothing here
          pos = surface->globalToLocal(Acts::GeometryContext(), {ahit.x(), ahit.y(), ahit.z()}, {0, 0, 0})
                    .value(); //, pos);

          Acts::BoundVector loc = Acts::BoundVector::Zero();
          loc[Acts::eBoundLoc0] = pos[0];
          loc[Acts::eBoundLoc0] = pos[1];
          debug() << "loc : (" << loc[0] << ", " << loc[1] << ")" << endmsg;

          // create source link at the end of the container
          auto it = source_links->emplace_hint(source_links->end(), *surface, 2, loc, cov);
          // ensure hits and links share the same order to prevent ugly surprises
          if (std::next(it) != source_links->end()) {
            error() << "The hit ordering broke. Run for your life." << endmsg;
            return StatusCode::FAILURE;
          }
        }
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(TrackingHitsSourceLinker)

} // namespace Jug::Reco
