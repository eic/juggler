#include <algorithm>
#include <cmath>

#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/RndmGenerators.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
#include "dd4pod/CalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitData.h"

namespace Jug {
  namespace Digi {

    /** Ecal Tungsten Sampling Calorimeter detector digitization.
     *
     *
     */
    class EcalTungstenSamplingDigi : public GaudiAlgorithm {
    public:
      Gaudi::Property<double>                      m_energyResolution{this, "energyResolution", 0.11}; // 11%sqrt(E)
      Rndm::Numbers                                m_gaussDist;
      DataHandle<dd4pod::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                        this};
      DataHandle<eic::RawCalorimeterHitCollection> m_outputHitCollection{"outputHitCollection",
                                                                         Gaudi::DataHandle::Writer, this};

      //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
      EcalTungstenSamplingDigi(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("inputHitCollection", m_inputHitCollection, "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
      }

      StatusCode initialize() override
      {
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;
        IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        StatusCode   sc      = m_gaussDist.initialize(randSvc, Rndm::Gauss(0.0, m_energyResolution.value()));
        if (!sc.isSuccess()) {
          return StatusCode::FAILURE;
        }
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override
      {
        // Note the energy is in units of GeV from dd4hep
        // input collections
        const dd4pod::CalorimeterHitCollection* simhits = m_inputHitCollection.get();
        // Create output collections
        auto                              rawhits          = m_outputHitCollection.createAndPut();
        eic::RawCalorimeterHitCollection* rawHitCollection = new eic::RawCalorimeterHitCollection();
        for (const auto& ahit : *simhits) {
          double                 res = m_gaussDist() / sqrt(ahit.energyDeposit());
          eic::RawCalorimeterHit rawhit(
              (long long)ahit.cellID(),
              std::llround(ahit.energyDeposit() * (1. + res) * 1.0e6), // convert to keV integer
              (double)ahit.truth().time * 1.0e6);
          rawhits->push_back(rawhit);
        }
        return StatusCode::SUCCESS;
      }
    };
    DECLARE_COMPONENT(EcalTungstenSamplingDigi)
  } // namespace Digi
} // namespace Jug
