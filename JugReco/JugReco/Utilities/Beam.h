#pragma once

#include "Math/Vector4D.h"
using ROOT::Math::PxPyPzEVector;

#include "edm4hep/MCParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"

namespace Jug::Reco {

namespace Beam {

  template<class collection>
  auto find_first_with_pdg(
      const collection& parts,
      const std::set<int32_t>& pdg) {
    collection c;
    for (const auto& p: parts) {
      if (pdg.count(p.getPDG()) > 0) {
        c.push_back(p);
        break;
      }
    }
    return c;
  }

  template<class collection>
  auto find_first_with_status_pdg(
      const collection& parts,
      const std::set<int32_t>& status,
      const std::set<int32_t>& pdg) {
    collection c;
    for (const auto& p: parts) {
      if (status.count(p.getGeneratorStatus()) > 0 && pdg.count(p.getPDG()) > 0) {
        c.push_back(p);
        break;
      }
    }
    return c;
  }

  auto find_first_beam_electron(const edm4hep::MCParticleCollection& mcparts) {
    return find_first_with_status_pdg(mcparts, {4}, {11});
  }

  auto find_first_beam_hadron(const edm4hep::MCParticleCollection& mcparts) {
    return find_first_with_status_pdg(mcparts, {4}, {2212, 2112});
  }

  auto find_first_scattered_electron(const edm4hep::MCParticleCollection& mcparts) {
    return find_first_with_status_pdg(mcparts, {1}, {11});
  }

  auto find_first_scattered_electron(const eicd::ReconstructedParticleCollection& rcparts) {
    return find_first_with_pdg(rcparts, {11});
  }

  PxPyPzEVector
  round_beam_four_momentum(
      const edm4hep::Vector3f& p_in,
      const float mass,
      const std::vector<float>& pz_set,
      const float crossing_angle = 0.0) {
    PxPyPzEVector p_out;
    for (const auto& pz : pz_set) {
      if (fabs(p_in.z / pz - 1) < 0.1) {
        p_out.SetPz(pz);
        break;
      }
    }
    p_out.SetPx(p_out.Pz() * sin(crossing_angle));
    p_out.SetPz(p_out.Pz() * cos(crossing_angle));
    p_out.SetE(std::hypot(p_out.Px(), p_out.Pz(), mass));
    return p_out;
  }

} // namespace Boost

} // namespace JugReco
