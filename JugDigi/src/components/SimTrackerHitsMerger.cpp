// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Whitney Armstrong, Chao Peng

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"

#include "JugBase/DataHandle.h"

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
    class SimTrackerHitsMerger : public GaudiAlgorithm {
    private:
      Gaudi::Property<std::vector<std::string>> m_inputSimTrackerHits{this, "inputSimTrackerHits", {},"Tracker hits to be aggregated"};
      DataHandle<edm4hep::SimTrackerHitCollection> m_outputSimTrackerHits{"outputSimTrackerHits", Gaudi::DataHandle::Writer, this};

      Gaudi::Property<std::vector<double>> m_inputProbabilities{this, "probabilities", {}};

      std::vector<DataHandle<edm4hep::SimTrackerHitCollection>*> m_hitCollections;

    public:
      SimTrackerHitsMerger(const std::string& name, ISvcLocator* svcLoc)
          : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("outputSimTrackerHits", m_outputSimTrackerHits, "output hits combined into single collection");
      }
      ~SimTrackerHitsMerger() {
        for (auto* col : m_hitCollections) {
          delete col;
        }
      }

      StatusCode initialize() override {
        if (GaudiAlgorithm::initialize().isFailure()) {
          return StatusCode::FAILURE;
        }
        for (auto colname : m_inputSimTrackerHits) {
          debug() << "initializing collection: " << colname  << endmsg;
          m_hitCollections.push_back(new DataHandle<edm4hep::SimTrackerHitCollection>{colname, Gaudi::DataHandle::Reader, this});
        }
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override
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
    DECLARE_COMPONENT(SimTrackerHitsMerger)

} // namespace Jug::Digi
