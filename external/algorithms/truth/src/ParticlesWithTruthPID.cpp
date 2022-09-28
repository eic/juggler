// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

#include <algorithms/truth/ParticlesWithTruthPID.h>

#include <algorithm>
#include <cmath>
#include <fmt/format.h>

#include <edm4eic/vector_utils.h>

namespace algorithms::truth {

void ParticlesWithTruthPID::init() {
  ; // do nothing
}
void ParticlesWithTruthPID::process(const ParticlesWithTruthPID::Input& input,
                                    const ParticlesWithTruthPID::Output& output) {
  const auto [mc_ptr, tracks_ptr] = input;
  auto [part_ptr, assoc_ptr]      = output;

  const auto& mc     = *mc_ptr;
  const auto& tracks = *tracks_ptr;
  auto& part         = *part_ptr;
  auto& assoc        = *assoc_ptr;

  const double sinPhiOver2Tolerance = sin(0.5 * m_phiTolerance);
  std::vector<bool> consumed(mc.size(), false);
  for (const auto& trk : tracks) {
    const auto mom =
        edm4eic::sphericalToVector(1.0 / std::abs(trk.getQOverP()), trk.getTheta(), trk.getPhi());
    const auto charge_rec = trk.getCharge();
    // utility variables for matching
    int best_match    = -1;
    double best_delta = std::numeric_limits<double>::max();
    for (size_t ip = 0; ip < mc.size(); ++ip) {
      const auto& mcpart = mc[ip];
      if (consumed[ip] || mcpart.getGeneratorStatus() > 1 || mcpart.getCharge() == 0 ||
          mcpart.getCharge() * charge_rec < 0) {
        if (aboveDebugThreshold()) {
          debug() << "ignoring non-primary/neutral/opposite charge particle" << endmsg;
        }
        continue;
      }
      const auto& p       = mcpart.getMomentum();
      const auto p_mag    = std::hypot(p.x, p.y, p.z);
      const auto p_phi    = std::atan2(p.y, p.x);
      const auto p_eta    = std::atanh(p.z / p_mag);
      const double dp_rel = std::abs((edm4eic::magnitude(mom) - p_mag) / p_mag);
      // check the tolerance for sin(dphi/2) to avoid the hemisphere problem and allow
      // for phi rollovers
      const double dsphi = std::abs(sin(0.5 * (edm4eic::angleAzimuthal(mom) - p_phi)));
      const double deta  = std::abs((edm4eic::eta(mom) - p_eta));

      if (dp_rel < m_pRelativeTolerance && deta < m_etaTolerance && dsphi < sinPhiOver2Tolerance) {
        const double delta = std::hypot(dp_rel / m_pRelativeTolerance, deta / m_etaTolerance,
                                        dsphi / sinPhiOver2Tolerance);
        if (delta < best_delta) {
          best_match = ip;
          best_delta = delta;
        }
      }
    }
    auto rec_part       = part.create();
    int32_t best_pid    = 0;
    auto referencePoint = rec_part.referencePoint();
    // float time          = 0;
    float mass = 0;
    if (best_match >= 0) {
      consumed[best_match] = true;
      const auto& mcpart   = mc[best_match];
      best_pid             = mcpart.getPDG();
      referencePoint       = {
          static_cast<float>(mcpart.getVertex().x), static_cast<float>(mcpart.getVertex().y),
          static_cast<float>(
              mcpart.getVertex().z)}; // @TODO: not sure if vertex/reference poitn makes sense here
      // time                 = mcpart.getTime();
      mass = mcpart.getMass();
    }
    rec_part.setType(static_cast<int16_t>(best_match >= 0 ? 0 : -1)); // @TODO: determine type codes
    rec_part.setEnergy(std::hypot(edm4eic::magnitude(mom), mass));
    rec_part.setMomentum(mom);
    rec_part.setReferencePoint(referencePoint);
    rec_part.setCharge(charge_rec);
    rec_part.setMass(mass);
    rec_part.setGoodnessOfPID(1); // perfect PID
    rec_part.setPDG(best_pid);
    // rec_part.covMatrix()  // @TODO: covariance matrix on 4-momentum
    // Also write MC <--> truth particle association if match was found
    if (best_match >= 0) {
      auto rec_assoc = assoc.create();
      rec_assoc.setRecID(rec_part.getObjectID().index);
      rec_assoc.setSimID(mc[best_match].getObjectID().index);
      rec_assoc.setWeight(1);
      rec_assoc.setRec(rec_part);
      // rec_assoc.setSim(mc[best_match]);
    }
    if (aboveDebugThreshold()) {
      if (best_match > 0) {
        const auto& mcpart = mc[best_match];
        debug() << fmt::format("Matched track with MC particle {}\n", best_match) << endmsg;
        debug() << fmt::format("  - Track: (mom: {}, theta: {}, phi: {}, charge: {})",
                               edm4eic::magnitude(mom), edm4eic::anglePolar(mom),
                               edm4eic::angleAzimuthal(mom), charge_rec)
                << endmsg;
        const auto& p      = mcpart.getMomentum();
        const auto p_mag   = edm4eic::magnitude(p);
        const auto p_phi   = edm4eic::angleAzimuthal(p);
        const auto p_theta = edm4eic::anglePolar(p);
        debug() << fmt::format(
                       "  - MC particle: (mom: {}, theta: {}, phi: {}, charge: {}, type: {}", p_mag,
                       p_theta, p_phi, mcpart.getCharge(), mcpart.getPDG())
                << endmsg;
      } else {
        debug() << fmt::format("Did not find a good match for track \n") << endmsg;
        debug() << fmt::format("  - Track: (mom: {}, theta: {}, phi: {}, charge: {})",
                               edm4eic::magnitude(mom), edm4eic::anglePolar(mom),
                               edm4eic::angleAzimuthal(mom), charge_rec)
                << endmsg;
      }
    }
  }
}

} // namespace algorithms::truth

