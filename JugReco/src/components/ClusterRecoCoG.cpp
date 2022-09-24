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

};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

