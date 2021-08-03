#include "JugTrack/GeometryContainers.hpp"
#include <ios>
#include "fmt/format.h"

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

#include "JugTrack/Index.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Measurement.hpp"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

  /** Sources Linker.
   *
   * \ingroup track
   */
  class TrackerSourcesLinker : public GaudiAlgorithm {
  public:

    Gaudi::Property<std::vector<std::string>>           u_inputHitCollections{this, "inputHitCollections", {}};
    std::vector<DataHandle<eic::TrackerHitCollection>*> u_hitCollections;
    DataHandle<IndexSourceLinkContainer> m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
    DataHandle<MeasurementContainer>     m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

  public:
    TrackerSourcesLinker(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("outputSourceLinks", m_outputSourceLinks, "");
      declareProperty("outputMeasurements", m_outputMeasurements, "");
    }

    ~TrackerSourcesLinker() {
      for (auto col : u_hitCollections) {
        if (col) { delete col; }
      }
    }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;

      for(const auto colname : u_inputHitCollections.value()) {
        auto new_col = new DataHandle<eic::TrackerHitCollection>{colname, Gaudi::DataHandle::Reader, this};
        u_hitCollections.push_back(new_col);
      }

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
      // Create output collections
      auto sourceLinks  = m_outputSourceLinks.createAndPut();
      auto measurements = m_outputMeasurements.createAndPut();

      int ihit = 0;
      for(const auto col : u_hitCollections) {
        for (const auto& ahit : *col->get()) {
          Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
          cov(0, 0) = ahit.covMatrix().covsym_xx * Acts::UnitConstants::mm * ahit.covMatrix().covsym_xx * Acts::UnitConstants::mm;
          cov(1, 1) = ahit.covMatrix().covsym_yy * Acts::UnitConstants::mm * ahit.covMatrix().covsym_yy * Acts::UnitConstants::mm;

          auto       vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
          auto       vol_id  = vol_ctx->identifier;
          const auto is      = m_geoSvc->surfaceMap().find(vol_id);
          if (is == m_geoSvc->surfaceMap().end()) {
            debug() << " vol_id (" << vol_id << ")  not found in m_surfaces." << endmsg;
            continue;
          }
          const Acts::Surface* surface = is->second;
          // transform global position into local coordinates
          Acts::Vector2 pos(0, 0);
          // geometry context contains nothing here
          pos = surface->globalToLocal(Acts::GeometryContext(), {ahit.position().x, ahit.position().y, ahit.position().z}, {0, 0, 0}).value();

          Acts::Vector2 loc     = Acts::Vector2::Zero();
          loc[Acts::eBoundLoc0] = pos[0]; //+ m_cfg.sigmaLoc0 * stdNormal(rng);
          loc[Acts::eBoundLoc1] = pos[1]; //+ m_cfg.sigmaLoc1 * stdNormal(rng);

          // the measurement container is unordered and the index under which the
          // measurement will be stored is known before adding it.
          Index           hitIdx = measurements->size();
          IndexSourceLink sourceLink(surface->geometryId(), hitIdx);
          auto            meas = Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);

          // add to output containers. since the input is already geometry-order,
          // new elements in geometry containers can just be appended at the end.
          sourceLinks->emplace_hint(sourceLinks->end(), std::move(sourceLink));
          measurements->emplace_back(std::move(meas));

          ihit++;
        }
      }
      debug() << fmt::format("{:d} hits linked from collections [{}]",
                             ihit, fmt::join(u_inputHitCollections.value(), ", "))
              << endmsg;
      return StatusCode::SUCCESS;
    }

   };
  DECLARE_COMPONENT(TrackerSourcesLinker)

} // namespace Jug::Reco
