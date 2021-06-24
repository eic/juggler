// Digitize the simulation outputs from Ecal Tungsten Sampling Calorimeter
#include <algorithm>
#include <cmath>

#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/RndmGenerators.h"

// FCCSW
#include "JugBase/DataHandle.h"

// Event Model related classes
#include "dd4pod/CalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitData.h"

using namespace Gaudi::Units;

namespace Jug {
  namespace Digi {

    /** Ecal Tungsten Sampling Calorimeter detector digitization.
     *
     *
     */
    class EcalTungstenSamplingDigi : public GaudiAlgorithm {
    public:
      Gaudi::Property<double>                      m_eRes{this, "energyResolution", 0.11}; // a%/sqrt(E/GeV)
      Gaudi::Property<std::vector<double>>         u_eRes{this, "energyResolutions", {}}; // a%/sqrt(E/GeV) + b% + c%/E
      Gaudi::Property<double>                      m_tRes{this, "timineResolution", 0.1*ns};
      Gaudi::Property<double>                      m_eUnit{this, "inputEnergyUnit", GeV};
      Gaudi::Property<double>                      m_tUnit{this, "inputTimeUnit", ns};
      Gaudi::Property<int>                         m_capADC{this, "capacityADC", 8096};
      Gaudi::Property<double>                      m_dyRangeADC{this, "dynamicRangeADC", 100*MeV};
      Gaudi::Property<int>                         m_pedMeanADC{this, "pedestalMean", 400};
      Gaudi::Property<double>                      m_pedSigmaADC{this, "pedestalSigma", 3.2};
      Rndm::Numbers                                m_normDist;
      DataHandle<dd4pod::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                        this};
      DataHandle<eic::RawCalorimeterHitCollection> m_outputHitCollection{"outputHitCollection",
                                                                         Gaudi::DataHandle::Writer, this};
      double res[3] = {0., 0., 0.};

      //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
      EcalTungstenSamplingDigi(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
      {
        declareProperty("inputHitCollection", m_inputHitCollection, "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
      }

      StatusCode initialize() override
      {
        warning() << "Deprecated algorithm for digi/reco, use Jug::Digi::CalorimeterHitDigi"
                     "and Jug::Reco::CalorimeterHitReco instead" << endmsg;
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;
        IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        StatusCode   sc      = m_normDist.initialize(randSvc, Rndm::Gauss(0.0, 1.0));
        if (!sc.isSuccess()) {
          return StatusCode::FAILURE;
        }
        // set energy resolution
        res[0] = m_eRes.value();
        for (size_t i = 0; i < u_eRes.size() && i < 3; ++i) {
            res[i] = u_eRes[i];
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
          double resval = std::pow(m_normDist()*res[0] / sqrt(ahit.energyDeposit()*m_eUnit/GeV), 2)
                        + std::pow(m_normDist()*res[1], 2)
                        + std::pow(m_normDist()*res[2] / (ahit.energyDeposit()*m_eUnit/GeV), 2);
          resval = std::sqrt(resval);
          double ped = m_pedMeanADC + m_normDist()*m_pedSigmaADC;
          long long adc = std::llround(ped + ahit.energyDeposit()*(1. + resval) * m_eUnit/m_dyRangeADC*m_capADC);
          eic::RawCalorimeterHit rawhit(
              (long long)ahit.cellID(),
              (adc > m_capADC ? m_capADC.value() : adc),
              (double)ahit.truth().time*m_tUnit/ns + m_normDist()*m_tRes/ns);
          rawhits->push_back(rawhit);
        }
        return StatusCode::SUCCESS;
      }
    };
    DECLARE_COMPONENT(EcalTungstenSamplingDigi)
  } // namespace Digi
} // namespace Jug
