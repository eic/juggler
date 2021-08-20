#include <algorithm>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
//#include "GaudiKernel/ToolHandle.h"
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "DD4hep/DD4hepUnits.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/Utilities/UniqueID.hpp"

// Event Model related classes
//#include "GaudiExamples/MyTrack.h"
#include "eicd/RawTrackerHitCollection.h"
#include "eicd/TrackerHitCollection.h"

#include "DD4hep/DD4hepUnits.h"

namespace Jug {
  namespace Reco {

    /** Tracker hit reconstruction.
     *
     * \ingroup reco
     */
    class TrackerHitReconstruction : public GaudiAlgorithm {
    private:
      // Unique identifier for this hit type, based on the algorithm name
      using HitClassificationType = decltype(eic::TrackerHitData().type);
      const HitClassificationType m_type;
    public:
      Gaudi::Property<double>                  m_timeResolution{this, "timeResolution", 10};
      Rndm::Numbers                            m_gaussDist;
      DataHandle<eic::RawTrackerHitCollection> m_inputHitCollection{
          "inputHitCollection", Gaudi::DataHandle::Reader, this};
      DataHandle<eic::TrackerHitCollection> m_outputHitCollection{"outputHitCollection",
                                                                  Gaudi::DataHandle::Writer, this};
      /// Pointer to the geometry service
      SmartIF<IGeoSvc> m_geoSvc;

    public:
      //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
      TrackerHitReconstruction(const std::string& name, ISvcLocator* svcLoc)
          : GaudiAlgorithm(name, svcLoc)
          , m_type{uniqueID<HitClassificationType>(name)}
      {
        declareProperty("inputHitCollection", m_inputHitCollection, "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
      }

      StatusCode initialize() override
      {
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;
        m_geoSvc = service("GeoSvc");
        if (!m_geoSvc) {
          error() << "Unable to locate Geometry Service. "
                  << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                  << endmsg;
          return StatusCode::FAILURE;
        }
        IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        StatusCode sc = m_gaussDist.initialize(randSvc, Rndm::Gauss(0.0, m_timeResolution.value()));
        if (!sc.isSuccess()) {
          return StatusCode::FAILURE;
        }
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override
      {
        constexpr auto mm = dd4hep::mm;
        // input collection
        const eic::RawTrackerHitCollection* rawhits = m_inputHitCollection.get();
        // Create output collections
        auto rec_hits = m_outputHitCollection.createAndPut();

        debug() << " raw hits size : " <<  std::size(*rawhits) << endmsg;
        for (const auto& ahit : *rawhits) {
          // debug() << "cell ID : " << ahit.cellID() << endmsg;
          auto pos = m_geoSvc->cellIDPositionConverter()->position(ahit.cellID());
          auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(ahit.cellID());
          // debug() << " dim size : " <<  std::size(dim) << endmsg;
          // for(const auto& s : dim ) {
          //  debug() << "a size : " <<  s << endmsg;
          //}
          // std::array<double,3> posarr; pos.GetCoordinates(posarr);
          // std::array<double,3> dimarr; dim.GetCoordinates(posarr);
          // eic::TrackerHit hit;
          eic::TrackerHit hit{(long long)ahit.cellID(),
                              ahit.ID(),
                              {pos.x() / mm, pos.y() / mm, pos.z() / mm, (float)ahit.time()/1000}, // mm, ns
                              {(dim[0]/mm)*(dim[0]/mm), (dim[1]/mm)*(dim[1]/mm), 0.0, 0.0},        // 
                              m_type,
                              (float)ahit.charge() / 1.0e6, // GeV
                              0.0f};
          rec_hits->push_back(hit);
        }
        return StatusCode::SUCCESS;
      }
    };
    DECLARE_COMPONENT(TrackerHitReconstruction)

  } // namespace Reco
} // namespace Jug

