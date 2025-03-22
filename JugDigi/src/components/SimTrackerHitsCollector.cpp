// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Whitney Armstrong, Chao Peng

// Gaudi
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include <k4FWCore/DataHandle.h>

// Event Model related classes
#include "edm4hep/SimTrackerHitCollection.h"

namespace Jug::Digi {

    /** Collect the tracking hits into a single collection.
     *
     * \param inputSimTrackerHits [in] vector of collection names
     * \param outputSimTrackerHits [out] hits combined into one collection.
     *
     * \ingroup digi
     */
    class SimTrackerHitsCollector : public Gaudi::Algorithm {
    private:
      Gaudi::Property<std::vector<std::string>> m_inputSimTrackerHits{this, "inputSimTrackerHits", {},"Tracker hits to be aggregated"};
      mutable DataHandle<edm4hep::SimTrackerHitCollection> m_outputSimTrackerHits{"outputSimTrackerHits", Gaudi::DataHandle::Writer, this};

      std::vector<DataHandle<edm4hep::SimTrackerHitCollection>*> m_hitCollections;

    public:
      SimTrackerHitsCollector(const std::string& name, ISvcLocator* svcLoc)
          : Gaudi::Algorithm(name, svcLoc)
      {
        declareProperty("outputSimTrackerHits", m_outputSimTrackerHits, "output hits combined into single collection");
      }
      ~SimTrackerHitsCollector() {
        for (auto* col : m_hitCollections) {
          delete col;
        }
      }

      StatusCode initialize() override {
        if (Gaudi::Algorithm::initialize().isFailure()) {
          return StatusCode::FAILURE;
        }
        for (auto colname : m_inputSimTrackerHits) {
          debug() << "initializing collection: " << colname  << endmsg;
          m_hitCollections.push_back(new DataHandle<edm4hep::SimTrackerHitCollection>{colname, Gaudi::DataHandle::Reader, this});
        }
        return StatusCode::SUCCESS;
      }

      StatusCode execute(const EventContext&) const override
      {
        auto* outputHits = m_outputSimTrackerHits.createAndPut();
        if (msgLevel(MSG::DEBUG)) {
          debug() << "execute collector" << endmsg;
        }
        for(const auto& hits: m_hitCollections) {
          const edm4hep::SimTrackerHitCollection* hitCol = hits->get();
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
    DECLARE_COMPONENT(SimTrackerHitsCollector)

} // namespace Jug::Digi
