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

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
JUGALGO_DEFINE_ALGORITHM(ClusterRecoCoG, algorithms::calorimetry::ClusterRecoCoG, Jug::Reco)
