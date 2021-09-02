#include <algorithm>
#include <cmath>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitData.h"
#include "dd4pod/CalorimeterHitCollection.h"

namespace Jug {
  namespace Digi {

    /** Hadronic Calorimeter Digitization.
     *
     * \f$ \sigma/E = a/\sqrt{E} \oplus b \f$
     *
     * \param a stochastic term (0.5 is 50%)
     * \param b constant term (0.05 is 5%)
     *
     * Resolution terms are added in quadrature. 
     * When digitizing they are assumed to be independent random variables and are sampled as such.
     *
     * \ingroup digi
     */
    class HadronicCalDigi : public GaudiAlgorithm, AlgorithmIDMixin<> {
    public:
      Gaudi::Property<double>                      m_energyResolution_a{this, "energyResolution_a", 0.5 /*50 percent*/};
      Gaudi::Property<double>                      m_energyResolution_b{this, "energyResolution_b", 0.05 /* 5 percent*/};
      Rndm::Numbers                                m_gaussDist_a;
      Rndm::Numbers                                m_gaussDist_b;
      DataHandle<dd4pod::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                        this};
      DataHandle<eic::RawCalorimeterHitCollection> m_outputHitCollection{"outputHitCollection",
                                                                         Gaudi::DataHandle::Writer, this};

    public:
      HadronicCalDigi(const std::string& name, ISvcLocator* svcLoc)
          : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
        declareProperty("inputHitCollection", m_inputHitCollection, "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
      }

      StatusCode initialize() override {
        warning() << "Deprecated algorithm for digi/reco, use Jug::Digi::CalorimeterHitDigi"
                     "and Jug::Reco::CalorimeterHitReco instead" << endmsg;
        IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        StatusCode   sc      = m_gaussDist_a.initialize(randSvc, Rndm::Gauss(0.0, m_energyResolution_a.value() ));
        if (!sc.isSuccess()) {
          return StatusCode::FAILURE;
        }
        sc      = m_gaussDist_b.initialize(randSvc, Rndm::Gauss(0.0, m_energyResolution_b.value()));
        if (!sc.isSuccess()) {
          return StatusCode::FAILURE;
        }
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;
        // f_counter = m_starting_value.value();
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override {
        // input collection
        const dd4pod::CalorimeterHitCollection* simhits = m_inputHitCollection.get();
        // Create output collections
        auto                              rawhits          = m_outputHitCollection.createAndPut();
        eic::RawCalorimeterHitCollection* rawHitCollection = new eic::RawCalorimeterHitCollection();
        int nhits = 0;
        for (const auto& ahit : *simhits) {
          // std::cout << ahit << "\n";
          double  sqrtE = std::sqrt(ahit.energyDeposit()) ;
          double aterm = m_gaussDist_a()*sqrtE;
          double bterm = ahit.energyDeposit()*m_gaussDist_b();
          // here 1000 is arbitrary scale factor
          eic::RawCalorimeterHit rawhit({nhits++, algorithmID()}, 
                                        (long long)ahit.cellID(),
                                        std::llround(ahit.energyDeposit() +aterm + bterm * 1000), 0);
          rawhits->push_back(rawhit);
        }
        return StatusCode::SUCCESS;
      }
    };
    DECLARE_COMPONENT(HadronicCalDigi)

  } // namespace Digi
} // namespace Jug

