// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Chao, Chao Peng, Whitney Armstrong

/*
 *  Reconstruct the cluster with Center of Gravity method
 *  Logarithmic weighting is used for mimicing energy deposit in transverse direction
 *
 *  Author: Sylvester Joosten, Chao Peng (ANL), 09/20/2022
 */

#include <JugAlgo/Algorithm.h>
#include <algorithms/calorimetry/ClusterRecoCoG.h>

#include "Gaudi/Property.h"

namespace Jug::Reco {

namespace {
  using AlgoBase = Jug::Algo::Algorithm<algorithms::calorimetry::ClusterRecoCoG>;
}

class ClusterRecoCoG : public AlgoBase {

public:
  ClusterRecoCoG(const std::string& name, ISvcLocator* svcLoc) : AlgoBase(name, svcLoc) {}

  virtual StatusCode configure() {
    setAlgoProp("samplingFraction", m_sampFrac.value());
    setAlgoProp("logWeightBase", m_logWeightBase.value());
    setAlgoProp("energyWeight", m_energyWeight.value());
    setAlgoProp("moduleDimZName", m_moduleDimZName.value());
    warning() << "DBG DBG type info for etabounds: " << typeid(m_enableEtaBounds.value()).name()
              << endmsg;
    setAlgoProp("enableEtaBounds", m_enableEtaBounds.value());
    return StatusCode::SUCCESS;
  }

private:
  Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};
  Gaudi::Property<double> m_logWeightBase{this, "logWeightBase", 3.6};
  Gaudi::Property<std::string> m_energyWeight{this, "energyWeight", "log"};
  Gaudi::Property<std::string> m_moduleDimZName{this, "moduleDimZName", ""};
  Gaudi::Property<bool> m_enableEtaBounds{this, "enableEtaBounds", false};
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

