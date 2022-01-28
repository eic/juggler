// Apply Birks Law to correct the energy deposit
// It uses the contributions member in dd4pod::CalorimeterHit, so simulation must enable storeCalorimeterContributions
//
// Author: Chao Peng
// Date: 09/29/2021

#include <algorithm>
#include <cmath>

#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"
#include "JugBase/IParticleSvc.h"

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
  class CalorimeterBirksCorr : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:

    // digitization settings
    Gaudi::Property<double> m_birksConstant{this, "birksConstant", 0.126*mm/MeV};

    DataHandle<dd4pod::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                      this};
    DataHandle<dd4pod::CalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                       this};

    SmartIF<IParticleSvc> m_pidSvc;
    // unitless conterpart of arguments
    double birksConstant;

    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    CalorimeterBirksCorr(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin{name, info()}
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }

      m_pidSvc = service("ParticleSvc");
      if (!m_pidSvc) {
        error() << "Unable to locate Particle Service! "
                << "Make sure you have ParticleSvc in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }

      // using juggler internal units (GeV, mm, radian, ns)
      birksConstant = m_birksConstant.value() / mm * GeV;

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      auto& ohits = *m_outputHitCollection.createAndPut();
      for (const auto& hit : *m_inputHitCollection.get()) {
        double lightyield = 0.;
        for (auto &truth : hit.contributions()) {
          const double charge = m_pidSvc->particle(truth.pdgID).charge;
          // some tolerance for precision
          if (std::abs(charge) > 1e-5) {
            lightyield += truth.deposit / (1. + truth.deposit / truth.length * birksConstant);
          }
        }
        auto ohit = ohits->create();
        ohit.cellID(hit.cellID());
        ohit.flag(hit.flag());
        ohit.g4ID(hit.g4ID());
        ohit.position(hit.position());
        ohit.truth(hit.truth());
        // replace energy deposit with Birks Law corrected value
        ohit.energyDeposit(lightyield);
      }
      return StatusCode::SUCCESS;
    }
  };
  DECLARE_COMPONENT(CalorimeterBirksCorr)

} // namespace Jug::Digi
