#include "GaudiAlg/GaudiAlgorithm.h"

#include "JugBase/DataHandle.h"

#include "edm4hep/MCParticleCollection.h"

class ReadTestConsumer : public GaudiAlgorithm {

public:
  ReadTestConsumer(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), m_genParticles("MCParticles", Gaudi::DataHandle::Reader, this) {
    declareProperty("genParticles", m_genParticles, "mc particles to read");
  }

  ~ReadTestConsumer(){};

  StatusCode initialize() {
    warning() << "This is a deprecated test algorithm" << endmsg;
    return GaudiAlgorithm::initialize(); }

  StatusCode execute() {
    // Read the input
    const edm4hep::MCParticleCollection* mcparticles = m_genParticles.get();

    // Does the reading work?
    debug() << mcparticles << endmsg;
    debug() << "MCParticle size: " << mcparticles->size() << endmsg;
    // counter for debug messages below
    //int cntr = 0;
    // Loop over all input particles
    //for (const auto& mcpart : *mcparticles) {
    //  if (10 > cntr++) {
    //    debug() << "time: " << mcpart.time << endmsg;
    //  }
    //}
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() { return GaudiAlgorithm::finalize(); }

private:
  /// Particles to read
  DataHandle<edm4hep::MCParticleCollection> m_genParticles;
};
DECLARE_COMPONENT(ReadTestConsumer)
