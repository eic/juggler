#include <algorithm>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Property.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
//#include "GaudiExamples/MyTrack.h"
#include "eicd/RawTrackerHitCollection.h"
#include "eicd/TrackerHitCollection.h"
#include "JugReco/SourceLinks.h"

namespace Jug {
  namespace Reco {
  
    /** Ultra-fast silicon detector digitization.
     *
     */
   class TrackerHitReconstruction : public GaudiAlgorithm {
   public:
    Gaudi::Property<double>      m_timeResolution{this, "timeResolution", 10};
    Rndm::Numbers m_gaussDist;
    DataHandle<eic::RawTrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::TrackerHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

   public:
    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    TrackerHitReconstruction(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc) {
          declareProperty("inputHitCollection", m_inputHitCollection,"");
          declareProperty("outputHitCollection", m_outputHitCollection, "");

        }

    StatusCode initialize() override {
      if (GaudiAlgorithm::initialize().isFailure())
        return StatusCode::FAILURE;
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
        << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
        return StatusCode::FAILURE;
      }
      IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
      StatusCode   sc      = m_gaussDist.initialize( randSvc, Rndm::Gauss(0.0, m_timeResolution.value()));
      if (!sc.isSuccess()) {
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      // input collection
      const eic::RawTrackerHitCollection* rawhits = m_inputHitCollection.get();
      // Create output collections
      auto rec_hits = m_outputHitCollection.createAndPut();
      for(const auto& ahit : *rawhits) {
        debug() << "cell ID : " << ahit.cellID() << endmsg;
        auto pos = m_geoSvc->cellIDPositionConverter()->position(ahit.cellID());
        auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(ahit.cellID());
        debug() << " dim size : " <<  std::size(dim) << endmsg;
        for(const auto& s : dim ) {
          debug() << "a size : " <<  s << endmsg;
        }
        //std::array<double,3> posarr; pos.GetCoordinates(posarr);
        //std::array<double,3> dimarr; dim.GetCoordinates(posarr);
        //eic::TrackerHit hit;
        eic::TrackerHit hit((long long)ahit.cellID(), (long long)ahit.cellID(), (long long)ahit.time(),
                            (float)ahit.charge() / 10000.0, (float)0.0, {{pos.x(), pos.y(),pos.z()}},{{dim[0],dim[1],0.0}});
        rec_hits->push_back(hit);
      }
      return StatusCode::SUCCESS;
    }

  };
  DECLARE_COMPONENT(TrackerHitReconstruction)

  } // namespace Examples
} // namespace Gaudi

