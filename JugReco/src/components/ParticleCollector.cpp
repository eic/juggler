// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"

#include "JugBase/DataHandle.h"

// Event Model related classes
#include "eicd/ReconstructedParticleCollection.h"

namespace Jug::Reco {

/** Collect the tracking hits into a single collection.
 *
 * \param inputParticles [in] vector of collection names
 * \param outputParticles [out] all particles into one collection.
 *
 * \ingroup reco
 */
class ParticleCollector : public GaudiAlgorithm {
public:
  Gaudi::Property<std::vector<std::string>> m_inputParticles{this, "inputParticles", {}, "Particles to be aggregated"};
  DataHandle<eicd::ReconstructedParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer,
                                                                     this};

  std::vector<DataHandle<eicd::ReconstructedParticleCollection>*> m_particleCollections;

public:
  ParticleCollector(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("outputParticles", m_outputParticles, "output particles combined into single collection");
  }
  ~ParticleCollector() {
    for (auto col : m_particleCollections) {
      if (col) {
        delete col;
      }
    }
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    for (auto colname : m_inputParticles) {
      debug() << "initializing collection: " << colname << endmsg;
      m_particleCollections.push_back(
          new DataHandle<eicd::ReconstructedParticleCollection>{colname, Gaudi::DataHandle::Reader, this});
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    auto output = m_outputParticles.createAndPut();
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
DECLARE_COMPONENT(ParticleCollector)

} // namespace Jug::Reco
