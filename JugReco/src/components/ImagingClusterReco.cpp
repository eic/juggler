// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Wouter Deconinck

#include <JugAlgo/Algorithm.h>
#include <edm4hep/utils/vector_utils.h> // FIXME remove when https://github.com/eic/EICrecon/pull/1482 merged
#include <EICrecon/algorithms/calorimetry/ImagingClusterReco.h>

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
JUGALGO_DEFINE_ALGORITHM(ImagingClusterReco, eicrecon::ImagingClusterReco, Jug::Reco)
