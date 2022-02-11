// Deprecated algorithm, as we can now properly store input collections in our output

#include <algorithm>
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/GaudiTool.h"
#include "Gaudi/Algorithm.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
//#include "GaudiExamples/MyTrack.h"
#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"
#include "dd4pod/TrackerHitCollection.h"


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
        declareProperty("inputCollection", m_inputHitCollection, "mcparticles");
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

      DataHandle<T_IN> m_inputHitCollection{"mcparticles", Gaudi::DataHandle::Reader, this};
      DataHandle<T_OUT> m_outputHitCollection{"genparticles", Gaudi::DataHandle::Writer, this};
    };

    using CalorimeterColCopier = InputCopier<dd4pod::CalorimeterHitCollection, dd4pod::CalorimeterHitCollection>;
    DECLARE_COMPONENT(CalorimeterColCopier)

    using TrackerColCopier = InputCopier<dd4pod::TrackerHitCollection, dd4pod::TrackerHitCollection>;
    DECLARE_COMPONENT(TrackerColCopier)

    using MCCopier = InputCopier<dd4pod::Geant4ParticleCollection, dd4pod::Geant4ParticleCollection>;
    DECLARE_COMPONENT(MCCopier)


    //class MCCopier : public GaudiAlgorithm {
    //public:
    //  //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    //  MCCopier(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    //  {
    //    declareProperty("inputCollection", m_inputHitCollection, "mcparticles");
    //    declareProperty("outputCollection", m_outputHitCollection, "genparticles");
    //  }
    //  StatusCode initialize() override
    //  {
    //    if (GaudiAlgorithm::initialize().isFailure())
    //      return StatusCode::FAILURE;
    //    // f_counter = m_starting_value.value();
    //    return StatusCode::SUCCESS;
    //  }
    //  StatusCode execute() override
    //  {
    //    // input collection
    //    const dd4pod::Geant4ParticleCollection* simhits = m_inputHitCollection.get();
    //    // Create output collections
    //    auto out_parts = m_outputHitCollection.createAndPut();
    //    // std::copy(std::begin(*simhits),std::end(*simhits),std::begin(*out_parts));
    //    for (const auto& ahit : *simhits) {
    //      // std::cout << ahit << "\n";
    //      // eic::RawCalorimeterHit rawhit((long long)ahit.cellID(),
    //      //                (long long)ahit.energyDeposit() * 100, 0);
    //      // dd4pod::Geant4Particle bpart = ahit;// dd4pod::Geant4ParticleCollection():
    //      out_parts->push_back(ahit.clone());
    //    }
    //    return StatusCode::SUCCESS;
    //  }

    //  DataHandle<dd4pod::Geant4ParticleCollection> m_inputHitCollection{"mcparticles", Gaudi::DataHandle::Reader, this};
    //  DataHandle<dd4pod::Geant4ParticleCollection> m_outputHitCollection{"genparticles", Gaudi::DataHandle::Writer,
    //                                                                     this};
    //};
    //DECLARE_COMPONENT(MCCopier)

  } // namespace Examples
} // namespace Gaudi

