// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"

#include "JugBase/DataHandle.h"

// Event Model related classes
#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

    /** Collect the tracking hits into a single collection.
     *
     * \ingroup reco
     */
    class TrackingHitsCollector : public GaudiAlgorithm {
    public:
      DataHandle<eic::TrackerHitCollection> m_trackerBarrelHits{"trackerBarrelHits", Gaudi::DataHandle::Reader, this};
      DataHandle<eic::TrackerHitCollection> m_trackerEndcapHits{"trackerEndcapHits", Gaudi::DataHandle::Reader, this};
      DataHandle<eic::TrackerHitCollection> m_vertexBarrelHits {"vertexBarrelHits" , Gaudi::DataHandle::Reader, this};
      DataHandle<eic::TrackerHitCollection> m_vertexEndcapHits {"vertexEndcapHits" , Gaudi::DataHandle::Reader, this};
      DataHandle<eic::TrackerHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};

    public:
      TrackingHitsCollector(const std::string& name, ISvcLocator* svcLoc)
          : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("trackerBarrelHits", m_trackerBarrelHits, "");
        declareProperty("trackerEndcapHits", m_trackerEndcapHits, "");
        declareProperty("vertexBarrelHits" , m_vertexBarrelHits , "");
        declareProperty("vertexEndcapHits" , m_vertexEndcapHits , "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
      }

      StatusCode initialize() override
      {
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override
      {
        const eic::TrackerHitCollection* trkBarrelHits = m_trackerBarrelHits.get();
        const eic::TrackerHitCollection* trkEndcapHits = m_trackerEndcapHits.get();
        const eic::TrackerHitCollection* vtxBarrelHits = m_vertexBarrelHits .get();
        const eic::TrackerHitCollection* vtxEndcapHits = m_vertexEndcapHits .get();
        auto outputHits = m_outputHitCollection.createAndPut();

        if(trkBarrelHits) {
          for (const auto& ahit : *trkBarrelHits) {
            outputHits->push_back(ahit.clone());
          }
        }
        if(trkEndcapHits) {
          for (const auto& ahit : *trkEndcapHits) {
            outputHits->push_back(ahit.clone());
          }
        }
        if(vtxBarrelHits) {
          for (const auto& ahit : *vtxBarrelHits) {
            outputHits->push_back(ahit.clone());
          }
        }
        if(vtxEndcapHits) {
          for (const auto& ahit : *vtxEndcapHits) {
            outputHits->push_back(ahit.clone());
          }
        }

        return StatusCode::SUCCESS;
      }
    };
    DECLARE_COMPONENT(TrackingHitsCollector)

} // namespace Jug::Reco
