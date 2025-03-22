// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck

// Gaudi
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include <k4FWCore/DataHandle.h>

// Event Model related classes
#include "edm4eic/TrackerHitCollection.h"

namespace Jug::Reco {

    /** Collect the tracking hits into a single collection.
     *
     * \param inputTrackingHits [in] vector of collection names
     * \param trackingHits [out] hits combined into one collection.
     *
     * \ingroup reco
     */
    class TrackingHitsCollector2 : public Gaudi::Algorithm {
    private:
      Gaudi::Property<std::vector<std::string>> m_inputTrackingHits{this, "inputTrackingHits", {},"Tracker hits to be aggregated"};
      mutable DataHandle<edm4eic::TrackerHitCollection> m_trackingHits{"trackingHits", Gaudi::DataHandle::Writer, this};

      mutable std::vector<DataHandle<const edm4eic::TrackerHitCollection>*> m_hitCollections;

    public:
      TrackingHitsCollector2(const std::string& name, ISvcLocator* svcLoc)
          : Gaudi::Algorithm(name, svcLoc)
      {
        declareProperty("trackingHits", m_trackingHits, "output hits combined into single collection");
      }
      ~TrackingHitsCollector2() {
        for (auto* col : m_hitCollections) {
          delete col;
        }
      }

      StatusCode initialize() override {
        if (Gaudi::Algorithm::initialize().isFailure()) {
          return StatusCode::FAILURE;
        }
        for (auto colname : m_inputTrackingHits) {
          debug() << "initializing collection: " << colname  << endmsg;
          m_hitCollections.push_back(new DataHandle<const edm4eic::TrackerHitCollection>{colname, Gaudi::DataHandle::Reader, this});
        }
        return StatusCode::SUCCESS;
      }

      StatusCode execute(const EventContext&) const override
      {
        auto* outputHits = m_trackingHits.createAndPut();
        if (msgLevel(MSG::DEBUG)) {
          debug() << "execute collector" << endmsg;
        }
        for(const auto& hits: m_hitCollections) {
          const edm4eic::TrackerHitCollection* hitCol = hits->get();
          if (msgLevel(MSG::DEBUG)) {
            debug() << "col n hits: " << hitCol->size() << endmsg;
          }
          for (const auto& ahit : *hitCol) {
            outputHits->push_back(ahit.clone());
          }
        }

        return StatusCode::SUCCESS;
      }
    };
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    DECLARE_COMPONENT_WITH_ID(TrackingHitsCollector2, "TrackingHitsCollector")
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    DECLARE_COMPONENT(TrackingHitsCollector2)

} // namespace Jug::Reco
