// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Whitney Armstrong

// Gaudi
#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"

#include <k4FWCore/DataHandle.h>

// Event Model related classes
#include "edm4eic/ReconstructedParticleCollection.h"

namespace Jug::Reco {

/** Collect the tracking hits into a single collection.
 *
 * \param inputParticles [in] vector of collection names
 * \param outputParticles [out] all particles into one collection.
 *
 * \ingroup reco
 */
class ParticleCollector : public Gaudi::Algorithm {
private:
  Gaudi::Property<std::vector<std::string>> m_inputParticles{this, "inputParticles", {}, "Particles to be aggregated"};
  mutable DataHandle<edm4eic::ReconstructedParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer,
                                                                     this};

  std::vector<DataHandle<edm4eic::ReconstructedParticleCollection>*> m_particleCollections;

public:
  ParticleCollector(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
    declareProperty("outputParticles", m_outputParticles, "output particles combined into single collection");
  }
  ~ParticleCollector() {
    for (auto* col : m_particleCollections) {
      delete col;
    }
  }

  StatusCode initialize() override {
    if (Gaudi::Algorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    for (auto colname : m_inputParticles) {
      debug() << "initializing collection: " << colname << endmsg;
      m_particleCollections.push_back(
          new DataHandle<edm4eic::ReconstructedParticleCollection>{colname, Gaudi::DataHandle::Reader, this});
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    auto* output = m_outputParticles.createAndPut();
    if (msgLevel(MSG::DEBUG)) {
      debug() << "execute collector" << endmsg;
    }
    for (const auto& list : m_particleCollections) {
      const auto& parts = *(list->get());
      if (msgLevel(MSG::DEBUG)) {
        debug() << "col n particles: " << parts.size() << endmsg;
      }
      for (const auto& part : parts) {
        output->push_back(part.clone());
      }
    }

    return StatusCode::SUCCESS;
  }
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ParticleCollector)

} // namespace Jug::Reco
