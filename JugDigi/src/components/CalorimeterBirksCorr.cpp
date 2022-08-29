// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Wouter Deconinck, Whitney Armstrong, Sylvester Joosten, Jihee Kim

// Apply Birks Law to correct the energy deposit
// It uses the contributions member in edm4hep::CalorimeterHit, so simulation must enable storeCalorimeterContributions
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
#include "JugBase/IParticleSvc.h"

// Event Model related classes
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4eic/RawCalorimeterHitCollection.h"
#include "edm4eic/RawCalorimeterHitData.h"


using namespace Gaudi::Units;

namespace Jug::Digi {

  /** Generic calorimeter hit digitiziation.
   *
   * \ingroup digi
   * \ingroup calorimetry
   */
  class CalorimeterBirksCorr : public GaudiAlgorithm {
  public:

    // digitization settings
    Gaudi::Property<double> m_birksConstant{this, "birksConstant", 0.126*mm/MeV};

    DataHandle<edm4hep::SimCalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                      this};
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                       this};

    SmartIF<IParticleSvc> m_pidSvc;
    // unitless conterpart of arguments
    double birksConstant{0};

    //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    CalorimeterBirksCorr(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
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
        error() << "Unable to locate Particle Service. "
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
        auto ohit = ohits->create(hit.getCellID(), hit.getEnergy(), hit.getPosition());
        double energy = 0.;
        for (const auto &c: hit.getContributions()) {
          ohit.addToContributions(c);
          const double charge = m_pidSvc->particle(c.getPDG()).charge;
          // some tolerance for precision
          if (std::abs(charge) > 1e-5) {
            // FIXME
            //energy += c.getEnergy() / (1. + c.getEnergy() / c.length * birksConstant);
            error() << "edm4hep::CaloHitContribution has no length field for Birks correction." << endmsg;
            return StatusCode::FAILURE;
          }
        }
        // replace energy deposit with Birks Law corrected value
        ohit.setEnergy(energy);
      }
      return StatusCode::SUCCESS;
    }
  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(CalorimeterBirksCorr)

} // namespace Jug::Digi
