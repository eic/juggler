// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten

#include <algorithms/algorithm.h>
#include <algorithms/random.h>

#include <DD4hep/DD4hepUnits.h>

// Event Model related classes
#include "edm4eic/RawTrackerHitCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimTrackerHitCollection.h"

namespace algorithms::tracking {

using SiliconTrackerDigiAlgorithm =
    Algorithm<Input<edm4hep::SimTrackerHitCollection>, Output<edm4eic::RawTrackerHitCollection>>;

/** Silicon detector digitization.
 *
 * \ingroup digi
 */
class SiliconTrackerDigi : public SiliconTrackerDigiAlgorithm {
private:
  Property<double> m_timeResolution{this, "timeResolution", 10,
                                    "Time resolution"}; // todo : add units
  Property<double> m_threshold{this, "threshold", 0. * dd4hep::keV, "Energy threshold"};

public:
  SiliconTrackerDigi(std::string_view name)
      : SiliconTrackerDigiAlgorithm{name,
                                    {"inputHitCollection"},
                                    {"outputHitColection"},
                                    "Create digitized tracker hits from simulated tracker hits "
                                    "based on time resolutions resolutions"} {}

  void init();
  void process(const Input&, const Output&);

  Generator m_rng = RandomSvc::instance().generator();
};

} // namespace algorithms::tracking
