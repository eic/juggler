#include <algorithm>
#include <cmath>
#include <fmt/format.h>

#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "eicd/TrackerHitCollection.h"
#include "eicd/ReconstructedParticleCollection.h"

namespace Jug::Reco {

class FarForwardParticles : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  DataHandle<eic::TrackerHitCollection> m_inputHitCollection{"FarForwardTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ReconstructedParticleCollection> m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer, this};

public:
  FarForwardParticles(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("inputCollection", m_inputHitCollection, "FarForwardTrackerHits");
    declareProperty("outputCollection", m_outputParticles, "ReconstructedParticles");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;

    return StatusCode::SUCCESS;
  }
  
  StatusCode execute() override {
    const auto& mc = *(m_inputHitCollection.get());
    auto& rc       = *(m_outputParticles.createAndPut());

    return StatusCode::SUCCESS;
  }

};

DECLARE_COMPONENT(FarForwardParticles)

} // namespace Jug::Reco


