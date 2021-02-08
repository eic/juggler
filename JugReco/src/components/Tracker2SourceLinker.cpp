#include "JugReco/GeometryContainers.hpp"
#include <ios>

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
#include "DD4hep/VolumeManager.h"
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

  class Tracker2SourceLinker : public GaudiAlgorithm {
  public:

    using HitCol = eic::TrackerHitCollection;

    DataHandle<HitCol> m_OTrackerBarrelHits{"TrackerBarrelHits", Gaudi::DataHandle::Reader, this};
    DataHandle<HitCol> m_OTrackerEndcapHits{"TrackerEndcapHits", Gaudi::DataHandle::Reader, this};
    DataHandle<HitCol> m_allTrackerHits{"allTrackerHits", Gaudi::DataHandle::Writer, this};

    std::vector<DataHandle<HitCol>*> m_trackerHitCollections;
    //Gaudi::Property<std::vector<std::string>> m_trackerHitCollections{this, "trackerHitCollections"};
    DataHandle<IndexSourceLinkContainer>     m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
    DataHandle<MeasurementContainer>          m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    /// Lookup container for hit surfaces that generate smeared hits
    std::unordered_map<uint64_t, const Acts::Surface*> m_surfaces;

  public:
    Tracker2SourceLinker(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("TrackerBarrelHits", m_OTrackerBarrelHits, "");
      declareProperty("TrackerEndcapHits", m_OTrackerEndcapHits, "");
      declareProperty("outputSourceLinks", m_outputSourceLinks, "");
      declareProperty("outputMeasurements", m_outputMeasurements, "");
      declareProperty("allTrackerHits", m_allTrackerHits, "");

    }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;

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
      auto sourceLinks  = m_outputSourceLinks.createAndPut();
      auto measurements = m_outputMeasurements.createAndPut();
      auto allHits      = m_allTrackerHits.createAndPut();
      // IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
      // IndexMultimap<Index> hitSimHitsMap;
      //int total_hits = 0;
      //for(const auto col : m_trackerHitCollections) {
      //  total_hits += col->get();
      //}
      //sourceLinks->reserve(hits->size());
      //measurements->reserve(hits->size());

      debug() << "   m_trackerHitCollections  size :  " << m_trackerHitCollections.size() << endmsg;

      int ihit = 0;
      for(const auto col : m_trackerHitCollections) {
        // input collection
        const eic::TrackerHitCollection* hits = col->get();
        // const eic::TrackerHitCollection* hits = get<eic::TrackerHitCollection>("/Event/"+ col);

        debug() << (*hits).size() << " hits " << endmsg;
        for (const auto& ahit : *hits) {

          allHits->push_back(ahit);

          Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
          cov(0, 0) = ahit.covsym_xx() * Acts::UnitConstants::mm * ahit.covsym_xx() * Acts::UnitConstants::mm;
          cov(1, 1) = ahit.covsym_yy() * Acts::UnitConstants::mm * ahit.covsym_yy() * Acts::UnitConstants::mm;

          auto       vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
          auto       vol_id  = vol_ctx->identifier;
          const auto is      = m_surfaces.find(vol_id);
          if (is == m_surfaces.end()) {
            debug() << " vol_id (" << vol_id << ")  not found in m_surfaces." << endmsg;
            continue;
          }
          const Acts::Surface* surface = is->second;
          debug() << " surface center : " << surface->center(Acts::GeometryContext()) << endmsg;
//=======
//
//          //allHits->push_back(ahit.clone());
//
//          Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
//          cov(0, 0) = ahit.covsym_xx() * Acts::UnitConstants::mm * ahit.covsym_xx() * Acts::UnitConstants::mm;
//          cov(1, 1) = ahit.covsym_yy() * Acts::UnitConstants::mm * ahit.covsym_yy() * Acts::UnitConstants::mm;
//
//          auto       vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
//          auto       vol_id  = vol_ctx->identifier;
//          const auto is      = m_surfaces.find(vol_id);
//
//          debug() << " cell ID : 0x" << std::hex << ahit.cellID() << std::dec << endmsg;
//          debug() << " vol  ID : 0x" << std::hex << vol_id << std::dec << endmsg;
//          if (is == m_surfaces.end()) {
//            auto detector = m_geoSvc->detector()->volumeManager().lookupDetector(vol_id);
//            debug() << " detector 0x" << std::hex << detector.id() <<std::dec << ", " << detector.path() << endmsg;
//            debug() << " vol_id   " << vol_id << " (0x" << std::hex << vol_id <<std::dec<< ")  not found in m_surfaces." << endmsg;
//            debug() << "hit " << ihit << endmsg;
//            continue;
//          }
//          const Acts::Surface* surface = is->second;
//          //debug() << " surface center : " << surface->center(Acts::GeometryContext()) << endmsg;
//>>>>>>> fitting
          // transform global position into local coordinates
          Acts::Vector2 pos(0, 0);
          // geometry context contains nothing here
          pos = surface->globalToLocal(Acts::GeometryContext(), {ahit.x(), ahit.y(), ahit.z()}, {0, 0, 0})
                    .value(); //, pos);

          Acts::Vector2 loc     = Acts::Vector2::Zero();
          loc[Acts::eBoundLoc0] = pos[0]; //+ m_cfg.sigmaLoc0 * stdNormal(rng);
          loc[Acts::eBoundLoc1] = pos[1]; //+ m_cfg.sigmaLoc1 * stdNormal(rng);
          // debug() << "loc : (" << loc[0] << ", " << loc[1] << ")" << endmsg;

          // local position
          // auto loc = {ahit.x(), ahit.y(), ahit.z()} - vol_ctx->volumePlacement().position()
          // debug() << " hit          : \n" <<  ahit << endmsg;
          // debug() << " cell ID : " << ahit.cellID() << endmsg;
          // debug() << " position : (" <<  ahit.position(0) << ", " <<  ahit.position(1) << ", "<<  ahit.position(2) <<
          // ") " << endmsg; debug() << " vol_id       : " <<  vol_id << endmsg; debug() << " placment pos : " <<
          // vol_ctx->volumePlacement().position() << endmsg;

          // the measurement container is unordered and the index under which the
          // measurement will be stored is known before adding it.
          Index           hitIdx = measurements->size();
          IndexSourceLink sourceLink(vol_id, ihit);
          auto            meas = Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);

          // add to output containers. since the input is already geometry-order,
          // new elements in geometry containers can just be appended at the end.
          sourceLinks->emplace_hint(sourceLinks->end(), std::move(sourceLink));
          measurements->emplace_back(std::move(meas));

          ihit++;
        }
      }
      return StatusCode::SUCCESS;
    }

   };
  DECLARE_COMPONENT(Tracker2SourceLinker)

} // namespace Jug::Reco
