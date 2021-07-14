#include "JugTrack/GeometryContainers.hpp"
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
#include "DD4hep/DD4hepUnits.h"


#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"

#include "JugTrack/Index.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/ProtoTrack.hpp"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

  /** Single Track source Linker and proto tracks.
   * This algorithm assumes only single track events.
   *
   */
  class SingleTrackSourceLinker : public GaudiAlgorithm {
  public:
    DataHandle<eic::TrackerHitCollection>    m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<IndexSourceLinkContainer>     m_outputSourceLinks{"outputSourceLinks", Gaudi::DataHandle::Writer, this};
    DataHandle<MeasurementContainer>         m_outputMeasurements{"outputMeasurements", Gaudi::DataHandle::Writer, this};
    DataHandle<ProtoTrackContainer>          m_outputProtoTracks{"outputProtoTracks", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    /// Lookup container for hit surfaces that generate smeared hits
    std::unordered_map<uint64_t, const Acts::Surface*> m_surfaces;

  public:
    SingleTrackSourceLinker(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputSourceLinks", m_outputSourceLinks, "");
      declareProperty("outputMeasurements", m_outputMeasurements, "");
      declareProperty("outputProtoTracks", m_outputProtoTracks, "");
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
      // input collection
      const eic::TrackerHitCollection* hits = m_inputHitCollection.get();
      // Create output collections
      auto sourceLinks = m_outputSourceLinks.createAndPut();
      auto measurements = m_outputMeasurements.createAndPut();
      auto protoTracks = m_outputProtoTracks.createAndPut();
      // IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
      // IndexMultimap<Index> hitSimHitsMap;
      sourceLinks->reserve(hits->size());
      measurements->reserve(hits->size());


      // assume single track --> one ProtoTrack
      ProtoTrack track;
      track.reserve((*hits).size());

      debug() << (*hits).size() << " hits " << endmsg;
      int ihit = 0;
      for(const auto& ahit : *hits) {

        track.emplace_back(ihit);

        auto       vol_ctx = m_geoSvc->cellIDPositionConverter()->findContext(ahit.cellID());
        auto       vol_id  = vol_ctx->identifier;
        const auto is      = m_geoSvc->surfaceMap().find(vol_id);
        if (is == m_surfaces.end()) {
          debug() << " vol_id (" <<  vol_id << ")  not found in m_surfaces!!!!" <<endmsg;
          continue;
        }
        const Acts::Surface* surface = is->second;

        // NOTE
        // Here is where the important hit and tracking geometry is connected.
        // A "Measurement" is constructed to for each hit which makes the connection to
        // the tracking surface and covariance matrix

        auto volman         = m_geoSvc->detector()->volumeManager();
        auto detelem        = volman.lookupDetElement(vol_id&0x3FFF);
        auto center         = surface->center(Acts::GeometryContext());
        auto gcenter        = surface->localToGlobal(Acts::GeometryContext(), {center(0), center(1)}, {0, 0, 0});
        dd4hep::Position global_position(ahit.x(), ahit.y(), ahit.z());
        // dd4hep::Position global_position(ahit.x() - center(0), ahit.y() - center(1), ahit.z() - center(2));
        auto local_position = detelem.nominal().worldToLocal(global_position);

        debug() << "===== Debugging hit =====" << endmsg;
        debug() << vol_ctx->volumePlacement().volIDs().str() << ", " << detelem.name() << endmsg;
        debug() << fmt::format("cell id: {:#064b}", ahit.cellID()) << endmsg;
        debug() << fmt::format("vol id:  {:#064b}", vol_id) << endmsg;
        debug() << " surface center : (" << center(0) << ", " << center(1) << ", " << center(2) << ")" << endmsg;
        debug() << " surface gcenter: (" << gcenter(0) << ", " << gcenter(1) << ", " << gcenter(2) << ")" << endmsg;
        debug() << " detelem center: " << detelem.nominal().localToWorld({0, 0, 0}) << endmsg;
        debug() << "global pos (" << global_position.x() << ", " << global_position.y() << ", "
                << global_position.z() << ")" << endmsg;
        debug() << "local  pos (" << local_position.x() << ", " << local_position.y() << ", "
                << local_position.z() << ")" << endmsg;
        auto acts_pos = surface->globalToLocal(Acts::GeometryContext(),
                {ahit.x()*Acts::UnitConstants::mm, ahit.y()*Acts::UnitConstants::mm, ahit.z()*Acts::UnitConstants::mm},
                {0, 0, 0}).value();//, pos);
        Acts::Vector2 pos(acts_pos[0], acts_pos[1]);
        Acts::Vector3 gpos = surface->localToGlobal(Acts::GeometryContext(), pos, {0, 0, 0});//, pos);
        debug() << " ACTS local position  : (" << acts_pos[0] << ", " << acts_pos[1] << ")"<< endmsg;
        debug() << " ACTS global position  : (" << gpos(0) << ", " << gpos(1) << ", " << gpos(2) << ")"<< endmsg;
        // construct the vector of measured parameters (2d position in this case)

        // construct the covariance matrix
        Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
        cov(0,0) = ahit.covsym_xx()*Acts::UnitConstants::mm*ahit.covsym_xx()*Acts::UnitConstants::mm;
        cov(1,1) = ahit.covsym_yy()*Acts::UnitConstants::mm*ahit.covsym_yy()*Acts::UnitConstants::mm;

        // Above we only consider the two position coordinates the comment below shows how to add time
        // which we will probably want to try later.
        //
        //Acts::SymMatrix3 cov;
        //cov << 0.05, 0., 0., 0., 0.05, 0., 0., 0., 900. * Acts::UnitConstants::ps * Acts::UnitConstants::ps;
        //Acts::Vector3 par(localX, localY, simHit.time());

        Index           hitIdx = ihit;
        IndexSourceLink sourceLink(vol_id, ihit);
        auto meas = Acts::makeMeasurement(sourceLink, pos, cov,
                                          Acts::eBoundLoc0, Acts::eBoundLoc1);//, Acts::eBoundTime);

        // TODO: check that local  to global is the same in dd4hep and acts.
        //
        //debug() << "ACTS surface center : " << surface->center(Acts::GeometryContext()) << endmsg;
        // transform global position into local coordinates
        // geometry context contains nothing here
        //Acts::Vector2 loc = Acts::Vector2::Zero();
        //loc[Acts::eBoundLoc0]     = pos[0] ;//+ m_cfg.sigmaLoc0 * stdNormal(rng);
        //loc[Acts::eBoundLoc1]     = pos[1] ;//+ m_cfg.sigmaLoc1 * stdNormal(rng);
        ////debug() << "loc : (" << loc[0] << ", " << loc[1] << ")" << endmsg;

        ////local position
        ////auto loc = {ahit.x(), ahit.y(), ahit.z()} - vol_ctx->volumePlacement().position()
        ////debug() << " hit          : \n" <<  ahit << endmsg;
        ////debug() << " cell ID : " << ahit.cellID() << endmsg;
        ////debug() << " position : (" <<  ahit.position(0) << ", " <<  ahit.position(1) << ", "<<  ahit.position(2) << ") " << endmsg;
        ////debug() << " vol_id       : " <<  vol_id << endmsg;
        ////debug() << " placment pos : " << vol_ctx->volumePlacement().position() << endmsg;

        //// the measurement container is unordered and the index under which the
        //// measurement will be stored is known before adding it.
        //auto meas = Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);

        // add to output containers. since the input is already geometry-order,
        // new elements in geometry containers can just be appended at the end.
        sourceLinks->emplace_hint(sourceLinks->end(), std::move(sourceLink));
        measurements->emplace_back(std::move(meas));

        ihit++;
      }
      // add proto track to the output collection
      protoTracks->emplace_back(std::move(track));
      return StatusCode::SUCCESS;
    }

   };
  DECLARE_COMPONENT(SingleTrackSourceLinker)

} // namespace Jug::reco

