// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @brief All track-related shared types.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"

#include <vector>

namespace Jug {

  /// (Reconstructed) track parameters e.g. close to the vertex.
  using TrackParameters = ::Acts::BoundTrackParameters;
  /// Container of reconstructed track states for multiple tracks.
  using TrackParametersContainer = std::vector<TrackParameters>;

} // namespace Jug
