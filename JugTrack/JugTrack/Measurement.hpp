#ifndef JugTrack_Measurement_HH
#define JugTrack_Measurement_HH

#include "Acts/EventData/Measurement.hpp"
#include "JugTrack/IndexSourceLink.hpp"

#include <cassert>
#include <vector>

namespace Jug {

  /// Variable measurement type that can contain all possible combinations.
  using Measurement = ::Acts::BoundVariantMeasurement<IndexSourceLink>;
  /// Container of measurements.
  ///
  /// In contrast to the source links, the measurements themself must not be
  /// orderable. The source links stored in the measurements are treated
  /// as opaque here and no ordering is enforced on the stored measurements.
  using MeasurementContainer = std::vector<Measurement>;

  /// Calibrator to convert an index source link to a measurement.
  class MeasurementCalibrator {
  public:
    /// Construct an invalid calibrator. Required to allow copying.
    MeasurementCalibrator() = default;
    /// Construct using a user-provided container to chose measurements from.
    MeasurementCalibrator(const MeasurementContainer& measurements) : m_measurements(&measurements) {}

    /// Find the measurement corresponding to the source link.
    ///
    /// @tparam parameters_t Track parameters type
    /// @param sourceLink Input source link
    /// @param parameters Input track parameters (unused)
    template <typename parameters_t>
    const Measurement& operator()(const IndexSourceLink& sourceLink, const parameters_t& /* parameters */) const
    {
      assert(m_measurements and "Undefined measurement container in DigitizedCalibrator");
      assert((sourceLink.index() < m_measurements->size()) and "Source link index is outside the container bounds");
      return (*m_measurements)[sourceLink.m_index];
    }

  private:
    // use pointer so the calibrator is copyable and default constructible.
    const MeasurementContainer* m_measurements = nullptr;
  };

} // namespace Jug

#endif
