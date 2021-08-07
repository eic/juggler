// A general digitization for CalorimeterHit from simulation
// 1. Smear energy deposit with a/sqrt(E/GeV) + b + c/E or a/sqrt(E/GeV) (relative value)
// 2. Digitize the energy with dynamic ADC range and add pedestal (mean +- sigma)
// 3. Time conversion with smearing resolution (absolute value)
//
// Author: Chao Peng
// Date: 06/02/2021

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

namespace Jug::Digi {

  /** Generic calorimeter hit digitiziation.
   *
   * \ingroup digi
   * \ingroup calorimetry
   */
  class CalorimeterHitDigi : public GaudiAlgorithm {
  public:
    // additional smearing resolutions
    Gaudi::Property<std::vector<double>> u_eRes{this, "energyResolutions", {}}; // a/sqrt(E/GeV) + b + c/(E/GeV)
    Gaudi::Property<double>              m_tRes{this, "timineResolution", 0.0 * ns};

    // digitization settings
    Gaudi::Property<int>    m_capADC{this, "capacityADC", 8096};
    Gaudi::Property<double> m_dyRangeADC{this, "dynamicRangeADC", 100 * MeV};
    Gaudi::Property<int>    m_pedMeanADC{this, "pedestalMean", 400};
    Gaudi::Property<double> m_pedSigmaADC{this, "pedestalSigma", 3.2};
    Rndm::Numbers           m_normDist;

    DataHandle<dd4pod::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                      this};
    DataHandle<eic::RawCalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                       this};
    // unitless counterparts of inputs
    double dyRangeADC, tRes, eRes[3] = {0., 0., 0.};

    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    CalorimeterHitDigi(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      // random number generator from service
      auto randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
      auto sc      = m_normDist.initialize(randSvc, Rndm::Gauss(0.0, 1.0));
      if (!sc.isSuccess()) {
        return StatusCode::FAILURE;
      }
      // set energy resolution numbers
      for (size_t i = 0; i < u_eRes.size() && i < 3; ++i) {
        eRes[i] = u_eRes[i];
      }

      // using juggler internal units (GeV, mm, radian, ns)
      dyRangeADC = m_dyRangeADC.value() / GeV;
      tRes       = m_tRes.value() / ns;

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto simhits = m_inputHitCollection.get();
      // Create output collections
      auto rawhits = m_outputHitCollection.createAndPut();
      int nhits = 0;
      for (const auto& ahit : *simhits) {
        // Note: juggler internal unit of energy is GeV
        double                 eResRel = std::sqrt(std::pow(m_normDist() * eRes[0] / sqrt(ahit.energyDeposit()), 2) +
                                   std::pow(m_normDist() * eRes[1], 2) +
                                   std::pow(m_normDist() * eRes[2] / (ahit.energyDeposit()), 2));
        double                 ped     = m_pedMeanADC + m_normDist() * m_pedSigmaADC;
        long long              adc = std::llround(ped + ahit.energyDeposit() * (1. + eResRel) / dyRangeADC * m_capADC);
        eic::RawCalorimeterHit rawhit(
            (long long)ahit.cellID(), 
            (adc > m_capADC.value() ? m_capADC.value() : adc),
            static_cast<int64_t>(1e6*(double)ahit.truth().time + m_normDist() * tRes), // @FIXME: this shouldn't be hardcoded, but should still be stored as an integer type
            nhits++);
        rawhits->push_back(rawhit);
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(CalorimeterHitDigi)

} // namespace Jug::Digi
