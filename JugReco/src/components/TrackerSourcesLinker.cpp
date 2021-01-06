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


#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"

#include "JugReco/SourceLinks.h"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

  class TrackerSourcesLinker : public GaudiAlgorithm {
  public:

    using HitCol = eic::TrackerHitCollection;

    DataHandle<HitCol> m_ITrackerBarrelHits{"ITrackerBarrelHits", Gaudi::DataHandle::Reader, this};
    DataHandle<HitCol> m_ITrackerEndcapHits{"ITrackerEndcapHits", Gaudi::DataHandle::Reader, this};
    DataHandle<HitCol> m_OTrackerBarrelHits{"OTrackerBarrelHits", Gaudi::DataHandle::Reader, this};
    DataHandle<HitCol> m_OTrackerEndcapHits{"OTrackerEndcapHits", Gaudi::DataHandle::Reader, this};

    std::vector<DataHandle<HitCol>*> m_trackerHitCollections;
    //Gaudi::Property<std::vector<std::string>> m_trackerHitCollections{this, "trackerHitCollections"};
    DataHandle<SourceLinkContainer>           m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    /// Lookup container for hit surfaces that generate smeared hits
    std::unordered_map<uint64_t, const Acts::Surface*> m_surfaces;

  public:
    TrackerSourcesLinker(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("ITrackerBarrelHits", m_ITrackerBarrelHits, "Inner tracker barrel hits");
      declareProperty("ITrackerEndcapHits", m_ITrackerEndcapHits, "");
      declareProperty("OTrackerBarrelHits", m_OTrackerBarrelHits, "");
      declareProperty("OTrackerEndcapHits", m_OTrackerEndcapHits, "");
      declareProperty("outputSourceLinks", m_outputSourceLinks, "");

    }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;

      m_trackerHitCollections.push_back(&m_ITrackerBarrelHits);
      m_trackerHitCollections.push_back(&m_ITrackerEndcapHits);
      m_trackerHitCollections.push_back(&m_OTrackerBarrelHits);
      m_trackerHitCollections.push_back(&m_OTrackerEndcapHits);

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
      // Create output collections
      auto source_links = m_outputSourceLinks.createAndPut();
      // setup local covariance
      // TODO add support for per volume/layer/module settings

      debug() << "   m_trackerHitCollections  size :  " << m_trackerHitCollections.size() << endmsg;

      for(const auto col : m_trackerHitCollections) {
        // input collection
        const eic::TrackerHitCollection* hits = col->get();
        //const eic::TrackerHitCollection* hits = get<eic::TrackerHitCollection>("/Event/"+ col);

        debug() << (*hits).size() << " hits " << endmsg;
        for(const auto& ahit : *hits) {

          Acts::BoundMatrix cov           = Acts::BoundMatrix::Zero();
          cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = ahit.covsym_xx()*Acts::UnitConstants::mm*ahit.covsym_xx()*Acts::UnitConstants::mm;
          cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = ahit.covsym_yy()*Acts::UnitConstants::mm*ahit.covsym_yy()*Acts::UnitConstants::mm;

          auto vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
          auto vol_id = vol_ctx->identifier;
          //debug() << " hit          : \n" <<  ahit << endmsg;
          //debug() << "cell ID : " << ahit.cellID() << endmsg;
          //debug() << " position : (" <<  ahit.position(0) << ", " <<  ahit.position(1) << ", "<<  ahit.position(2) << ") " << endmsg;
          debug() << " vol_id       : " <<  vol_id << endmsg;
          debug() << " placment pos : " << vol_ctx->volumePlacement().position() << endmsg;

          const auto is = m_surfaces.find(vol_id);
          if (is == m_surfaces.end()) {
            debug() << " vol_id (" <<  vol_id << ")  not found in m_surfaces." <<endmsg;
            continue;
          }
          const Acts::Surface* surface = is->second;

          debug() << " surface center : " << surface->center(Acts::GeometryContext()) << endmsg;
          // transform global position into local coordinates
          Acts::Vector2D pos(0, 0);
          // geometry context contains nothing here
          pos = surface->globalToLocal(Acts::GeometryContext(), {ahit.x(), ahit.y(), ahit.z()}, {0, 0, 0}).value();//, pos);

          //// smear truth to create local measurement
          Acts::BoundVector loc = Acts::BoundVector::Zero();
          loc[Acts::eBoundLoc0]     = pos[0] ;//+ m_cfg.sigmaLoc0 * stdNormal(rng);
          loc[Acts::eBoundLoc0]     = pos[1] ;//+ m_cfg.sigmaLoc1 * stdNormal(rng);

          debug() << "loc : (" << loc[0] << ", " << loc[1] << ")" << endmsg;

          // create source link at the end of the container
          //auto it = source_links->emplace_hint(source_links->end(), *surface, hit, 2, loc, cov);
          auto it = source_links->emplace_hint(source_links->end(), *surface,  2, loc, cov);
          // ensure hits and links share the same order to prevent ugly surprises
          if (std::next(it) != source_links->end()) {
            error() << "The hit ordering broke. Run for your life." << endmsg;
            return StatusCode::FAILURE;
          }

          //std::array<double,3> posarr; pos.GetCoordinates(posarr);
          //std::array<double,3> dimarr; dim.GetCoordinates(posarr);
          //eic::TrackerHit hit;
          //eic::TrackerHit hit((long long)ahit.cellID0(), (long long)ahit.cellID(), (long long)ahit.time(),
          //                    (float)ahit.charge() / 10000.0, (float)0.0, {{pos.x(), pos.y(),pos.z()}},{{dim[0],dim[1],0.0}});
          //rec_hits->push_back(hit);
        }
      }
      return StatusCode::SUCCESS;
    }

   };
  DECLARE_COMPONENT(TrackerSourcesLinker)

} // namespace Jug::Reco
