#ifndef JUG_RECO_SourceLinks_HH
#define JUG_RECO_SourceLinks_HH

#include "Acts/EventData/Measurement.hpp"
#include "JugReco/GeometryContainers.hpp"

#include <stdexcept>
#include <string>
#include "dd4pod/Geant4Particle.h"

#include "eicd/TrackerHitCollection.h"


namespace Jug {

/** Source Link for simulation in the acts-framework.
 *
 * https://github.com/acts-project/acts/blob/master/Core/include/Acts/EventData/SourceLinkConcept.hpp
 * The source link stores the measuremts, surface, and the associated simulated
 * truth hit.
 *
 * @todo Allow multiple truth hits e.g. for merged hits.
 *
 */
class SourceLink {

 private:
  Acts::BoundVector m_values;
  Acts::BoundMatrix m_cov;
  size_t m_dim = 2;
  // store geo id copy to avoid indirection via truth hit
  int32_t m_index;
  Acts::GeometryIdentifier m_geometryId;
  // need to store pointers to make the object copyable
  const Acts::Surface* m_surface;
  //const ActsFatras::Hit* m_truthHit;
  const eic::TrackerHit* m_Hit ;

 public:
  SourceLink(const Acts::Surface& surface, //const ActsFatras::Hit& truthHit,
                size_t dim, int32_t index, Acts::BoundVector values, Acts::BoundMatrix cov)
      : m_values(values),
        m_cov(cov),
        m_dim(dim),
        m_index(index),
        m_geometryId(0),//truthHit.geometryId()),
        m_surface(&surface){}
        //m_truthHit(&truthHit) {}
  /// Must be default_constructible to satisfy SourceLinkConcept.
  SourceLink() = default;
  SourceLink(SourceLink&&) = default;
  SourceLink(const SourceLink&) = default;
  SourceLink& operator=(SourceLink&&) = default;
  SourceLink& operator=(const SourceLink&) = default;

  constexpr Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
  constexpr const Acts::Surface& referenceSurface() const { return *m_surface; }
  //constexpr const ActsFatras::Hit& truthHit() const { return *m_truthHit; }

  Acts::FittableMeasurement<SourceLink> operator*() const {
    if (m_dim == 0) {
      throw std::runtime_error("Cannot create dim 0 measurement");
    } else if (m_dim == 1) {
      return Acts::Measurement<SourceLink, Acts::BoundIndices,
                               Acts::eBoundLoc0>{
          m_surface->getSharedPtr(), *this, m_cov.topLeftCorner<1, 1>(),
          m_values[0]};
    } else if (m_dim == 2) {
      return Acts::Measurement<SourceLink, Acts::BoundIndices,
                               Acts::eBoundLoc0, Acts::eBoundLoc1>{
          m_surface->getSharedPtr(), *this, m_cov.topLeftCorner<2, 2>(),
          m_values[0], m_values[1]};
    } else {
      throw std::runtime_error("Dim " + std::to_string(m_dim) +
                               " currently not supported.");
    }
  }

  friend constexpr bool operator==(const SourceLink& lhs,
                                   const SourceLink& rhs) {

    return (lhs.geometryId() == rhs.geometryId()) && (lhs.m_index == rhs.m_index);
    //lhs.m_truthHit == rhs.m_truthHit;
  }
  friend constexpr bool operator!=(const SourceLink& lhs,
                                   const SourceLink& rhs) {
    return not(lhs == rhs);
  }
};

/// Store source links ordered by geometry identifier.
using SourceLinkContainer = GeometryIdMultiset<SourceLink>;
} // namespace Jug

#endif

