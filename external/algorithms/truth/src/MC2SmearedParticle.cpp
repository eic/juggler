// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Whitney Armstrong, Wouter Deconinck

#include <algorithms/truth/MC2SmearedParticle.h>

#include <cmath>
#include <edm4eic/vector_utils.h>

namespace algorithms::truth {

void MC2SmearedParticle::init() {
  ; // do nothing
}

void MC2SmearedParticle::process(const MC2SmearedParticle::Input& input,
                                 const MC2SmearedParticle::Output& output) const {
  const auto [parts] = input;
  auto [out_parts]   = output;

  for (const auto& p : *parts) {
    if (p.getGeneratorStatus() > 1) {
      if (aboveDebugThreshold()) {
        debug() << "ignoring particle with generatorStatus = " << p.getGeneratorStatus() << endmsg;
      }
      continue;
    }

    // for now just use total momentum smearing as this is the largest effect,
    // ideally we should also smear the angles but this should be good enough
    // for now.
    const auto pvec     = p.getMomentum();
    const auto pgen     = std::hypot(pvec.x, pvec.y, pvec.z);
    const auto momentum = pgen * m_rng->gaussian<double>(0., m_smearing);
    // make sure we keep energy consistent
    using MomType = decltype(edm4eic::ReconstructedParticle().getMomentum().x);
    const MomType energy =
        std::sqrt(p.getEnergy() * p.getEnergy() - pgen * pgen + momentum * momentum);
    const MomType px = p.getMomentum().x * momentum / pgen;
    const MomType py = p.getMomentum().y * momentum / pgen;
    const MomType pz = p.getMomentum().z * momentum / pgen;

    const MomType dpx = m_smearing * px;
    const MomType dpy = m_smearing * py;
    const MomType dpz = m_smearing * pz;
    const MomType dE  = m_smearing * energy;
    // ignore covariance for now
    // @TODO: vertex smearing
    const MomType vx = p.getVertex().x;
    const MomType vy = p.getVertex().y;
    const MomType vz = p.getVertex().z;

    auto rec_part = out_parts->create();
    rec_part.setType(-1); // @TODO: determine type codes
    rec_part.setEnergy(energy);
    rec_part.setMomentum({px, py, pz});
    rec_part.setReferencePoint({vx, vy, vz}); // @FIXME: probably not what we want?
    rec_part.setCharge(p.getCharge());
    rec_part.setMass(p.getMass());
    rec_part.setGoodnessOfPID(1); // Perfect PID
    rec_part.setCovMatrix({dpx, dpy, dpz, dE});
    rec_part.setPDG(p.getPDG());
  }
}

} // namespace algorithms::truth

