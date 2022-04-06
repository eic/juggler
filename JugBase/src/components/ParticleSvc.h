// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten
#ifndef PARTICLESVC_H
#define PARTICLESVC_H

#include <map>

#include "GaudiKernel/Service.h"

#include "JugBase/IParticleSvc.h"

/** Simple particle service.
 *
 *  This meant to provide basic particle information for reconstruction purposes.
 *  If particle data is needed, be sure to grab everything you can in an initialization
 *  step. Currently the  returned Particle/ParticleMap are by value.
 */
class ParticleSvc : public extends<Service, IParticleSvc> {
public:
  using Particle    = Jug::Base::ParticleData;
  using ParticleMap = std::map<int, Particle>;

  const ParticleMap m_particleMap;

public:
  ParticleSvc(const std::string& name, ISvcLocator* svc);

  virtual ~ParticleSvc();

  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final { return StatusCode::SUCCESS; }

  virtual ParticleMap particleMap() const { return m_particleMap; }
  virtual Particle particle(int pdg) const {
    if (m_particleMap.count(pdg) == 0) {
      // error
      return m_particleMap.at(0);
    }
    return m_particleMap.at(pdg);
  }
};

#endif
