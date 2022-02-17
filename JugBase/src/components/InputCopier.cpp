// Deprecated algorithm, as we can now properly store input collections in our output

#include <algorithm>
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/GaudiTool.h"
#include "Gaudi/Algorithm.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"


namespace Jug {
  namespace Base {

    /** Need to fix a bug.
     *
     * Details found here:
     * https://github.com/AIDASoft/podio/issues/103
     */
    template<typename T_IN, typename T_OUT>
    class InputCopier : public GaudiAlgorithm {
    public:
      InputCopier(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("inputCollection", m_inputHitCollection, "MCParticles");
        declareProperty("outputCollection", m_outputHitCollection, "genparticles");
      }
      StatusCode initialize() override
      {
        if (GaudiAlgorithm::initialize().isFailure()) {
          return StatusCode::FAILURE;
        }
        warning() << "DEPRECATED ALGORITHM, no need to use this anymore, we can do a proper straight passthrough from input to output." << endmsg;
        warning() << "1) Remove the calls to InputCopier from your options file." << endmsg;
        warning() << "2) Add 'keep mcparticles' to the PodioOutput.outputCommands." << endmsg;
        warning() << "3) Update your analysis code to use 'mcparticles' directly." << endmsg;
        return StatusCode::SUCCESS;
      }
      StatusCode execute() override
      {
        // input collection
        const T_IN* simhits = m_inputHitCollection.get();
        // output collection
        auto out_parts = m_outputHitCollection.createAndPut();
        for (const auto& ahit : *simhits) {
          out_parts->push_back(ahit.clone());
        }
        return StatusCode::SUCCESS;
      }

      DataHandle<T_IN> m_inputHitCollection{"MCParticles", Gaudi::DataHandle::Reader, this};
      DataHandle<T_OUT> m_outputHitCollection{"genparticles", Gaudi::DataHandle::Writer, this};
    };

    using CalorimeterColCopier = InputCopier<edm4hep::CalorimeterHitCollection, edm4hep::CalorimeterHitCollection>;
    DECLARE_COMPONENT(CalorimeterColCopier)

    using TrackerColCopier = InputCopier<edm4hep::TrackerHitCollection, edm4hep::TrackerHitCollection>;
    DECLARE_COMPONENT(TrackerColCopier)

    using MCCopier = InputCopier<edm4hep::MCParticleCollection, edm4hep::MCParticleCollection>;
    DECLARE_COMPONENT(MCCopier)

  } // namespace Examples
} // namespace Gaudi

