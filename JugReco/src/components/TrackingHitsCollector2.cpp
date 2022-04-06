// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck

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
    private:
      Gaudi::Property<std::vector<std::string>> m_inputTrackingHits{this, "inputTrackingHits", {},"Tracker hits to be aggregated"};
      DataHandle<eicd::TrackerHitCollection> m_trackingHits{"trackingHits", Gaudi::DataHandle::Writer, this};

      std::vector<DataHandle<eicd::TrackerHitCollection>*> m_hitCollections;

    public:
      TrackingHitsCollector2(const std::string& name, ISvcLocator* svcLoc)
          : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("trackingHits", m_trackingHits, "output hits combined into single collection");
      }
      ~TrackingHitsCollector2() {
        for (auto* col : m_hitCollections) {
          delete col;
        }
      }

      StatusCode initialize() override {
        if (GaudiAlgorithm::initialize().isFailure()) {
          return StatusCode::FAILURE;
        }
        for (auto colname : m_inputTrackingHits) {
          debug() << "initializing collection: " << colname  << endmsg;
          m_hitCollections.push_back(new DataHandle<eicd::TrackerHitCollection>{colname, Gaudi::DataHandle::Reader, this});
        }
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override
      {
        auto* outputHits = m_trackingHits.createAndPut();
        if (msgLevel(MSG::DEBUG)) {
          debug() << "execute collector" << endmsg;
        }
        for(const auto& hits: m_hitCollections) {
          const eicd::TrackerHitCollection* hitCol = hits->get();
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
    DECLARE_COMPONENT(TrackingHitsCollector2)

} // namespace Jug::Reco
