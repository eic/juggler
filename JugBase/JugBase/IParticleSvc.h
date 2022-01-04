#ifndef IParticleSvc_H
#define IParticleSvc_H

#include <GaudiKernel/IService.h>
#include <unordered_map>

namespace Jug::Base {

  /** Simple particle data.
   *
   */
  struct ParticleData {
    int         pdgCode;
    int         charge;
    double      mass; //std::string name;
  };
} // namespace Jug::Base

/** Particle interface.
 *
 * \ingroup base
 */
class GAUDI_API IParticleSvc : virtual public IService {
public:
  using Particle    = Jug::Base::ParticleData;
  using ParticleMap = std::map<int, Particle>;

public:
  /// InterfaceID
  DeclareInterfaceID(IParticleSvc, 1, 0);
  virtual ~IParticleSvc() {}

  virtual ParticleMap particleMap() const = 0;
  virtual Particle    particle(int pdg) const = 0;
};

#endif  // IParticleSvc_H
