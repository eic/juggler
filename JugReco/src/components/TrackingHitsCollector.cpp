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
      DataHandle<eicd::TrackerHitCollection> m_trackerBarrelHits{"trackerBarrelHits", Gaudi::DataHandle::Reader, this};
      DataHandle<eicd::TrackerHitCollection> m_trackerEndcapHits{"trackerEndcapHits", Gaudi::DataHandle::Reader, this};
      DataHandle<eicd::TrackerHitCollection> m_vertexBarrelHits {"vertexBarrelHits" , Gaudi::DataHandle::Reader, this};
      DataHandle<eicd::TrackerHitCollection> m_vertexEndcapHits {"vertexEndcapHits" , Gaudi::DataHandle::Reader, this};
      DataHandle<eicd::TrackerHitCollection> m_gemEndcapHits {"gemEndcapHits" , Gaudi::DataHandle::Reader, this};
      DataHandle<eicd::TrackerHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};

    public:
      TrackingHitsCollector(const std::string& name, ISvcLocator* svcLoc)
          : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("trackerBarrelHits", m_trackerBarrelHits, "");
        declareProperty("trackerEndcapHits", m_trackerEndcapHits, "");
        declareProperty("vertexBarrelHits" , m_vertexBarrelHits , "");
        declareProperty("vertexEndcapHits" , m_vertexEndcapHits , "");
        declareProperty("gemEndcapHits" , m_gemEndcapHits , "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
      }

      StatusCode initialize() override
      {
        warning() << "DEPRECATED, use TrackingHitsCollector2 instead" << endmsg;
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override
      {
        const eicd::TrackerHitCollection* trkBarrelHits = m_trackerBarrelHits.get();
        const eicd::TrackerHitCollection* trkEndcapHits = m_trackerEndcapHits.get();
        const eicd::TrackerHitCollection* vtxBarrelHits = m_vertexBarrelHits .get();
        const eicd::TrackerHitCollection* vtxEndcapHits = m_vertexEndcapHits .get();
        const eicd::TrackerHitCollection* gemEndcapHits = m_gemEndcapHits .get();
        auto outputHits = m_outputHitCollection.createAndPut();

        for (const auto* hits : {trkBarrelHits, trkEndcapHits, vtxBarrelHits, vtxEndcapHits, gemEndcapHits}) {
          if (hits) {
            for (const auto& ahit : *hits) {
              auto new_hit = ahit.clone();
              outputHits->push_back(ahit.clone());
            }
          }
        }

        return StatusCode::SUCCESS;
      }
    };
    DECLARE_COMPONENT(TrackingHitsCollector)

} // namespace Jug::Reco
