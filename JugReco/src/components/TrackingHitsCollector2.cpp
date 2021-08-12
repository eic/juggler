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
     * \param inputTrackingHits [in] vector of collection names
     * \param trackingHits [out] hits combined into one collection.
     *
     * \ingroup reco
     */
    class TrackingHitsCollector2 : public GaudiAlgorithm {
    public:
      Gaudi::Property<std::vector<std::string>> m_inputTrackingHits{this, "inputTrackingHits", {}};
      DataHandle<eic::TrackerHitCollection> m_trackingHits{"trackingHits", Gaudi::DataHandle::Writer, this};

      std::vector<std::unique_ptr<DataHandle<eic::TrackerHitCollection>>> m_hitCollections;

    public:
      TrackingHitsCollector2(const std::string& name, ISvcLocator* svcLoc)
          : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("inputTrackingHits" , m_inputTrackingHits , "vector of collection names");
        declareProperty("trackingHits", m_trackingHits, "output hits combined into single collection");
      }

      StatusCode initialize() override {
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;
        for (const auto colname : m_inputTrackingHits.value()) {
          m_hitCollections.push_back(std::make_unique<DataHandle<eic::TrackerHitCollection>>(colname, Gaudi::DataHandle::Reader, this));
        }
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override
      {
        auto outputHits = m_trackingHits.createAndPut();

        for(const auto& hits: m_hitCollections) {
          const eic::TrackerHitCollection* hitCol = hits->get();
          for (const auto& ahit : *hitCol) {
            outputHits->push_back(ahit.clone());
          }
        }

        return StatusCode::SUCCESS;
      }
    };
    DECLARE_COMPONENT(TrackingHitsCollector2)

} // namespace Jug::Reco
