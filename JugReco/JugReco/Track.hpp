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

#include "JugReco/SimMultiTrajectory.hpp"
//#include "ACTFW/EventData/SimSourceLink.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "JugReco/IndexSourceLink.hpp"

#include <vector>

namespace Jug {

  /// (Reconstructed) track parameters e.g. close to the vertex.
  using TrackParameters = ::Acts::BoundTrackParameters;
  /// Container of reconstructed track states for multiple tracks.
  using TrackParametersContainer = std::vector<TrackParameters>;

  /// MultiTrajectory definition
  using Trajectory = Acts::MultiTrajectory<IndexSourceLink>;

  /// Container for the truth fitting/finding track(s)
  using TrajectoryContainer = std::vector<SimMultiTrajectory>;

} // namespace Jug
