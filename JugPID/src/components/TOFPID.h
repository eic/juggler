
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/IParticleSvc.h"
#include "JugBase/UniqueID.h"
#include "JugBase/BField/DD4hepBField.h"
#include "JugTrack/Trajectories.hpp"

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/TrackerHitCollection.h"

#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/TofParticleIDCollection.h"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
//#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"

using BoundTrackParamPtr = std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;

#ifndef _TOF_PID_
#define _TOF_PID_

namespace Jug::PID {
  /** RICH IRT
   *
   * \ingroup reco
   */
  class TOFPID : public GaudiAlgorithm, AlgorithmIDMixin<> {

  public:
    // constructor
    TOFPID(const std::string& name, ISvcLocator* svcLoc);

    // INITIALIZE .....................................................
    StatusCode initialize( void ) override;

    // EXECUTE .....................................................
    StatusCode execute( void ) override;

    // FINALIZE .....................................................
    StatusCode finalize( void ) override;

    // Input collections;
    std::unique_ptr<DataHandle<dd4pod::TrackerHitCollection>> m_inputHitCollection;
    DataHandle<dd4pod::Geant4ParticleCollection> m_inputMCParticles {
      "inputMCParticles", 
	Gaudi::DataHandle::Reader, 
	this
	};
    DataHandle<eic::ReconstructedParticleCollection> m_inputRecoParticles {
      "inputRecoParticles", 
	Gaudi::DataHandle::Reader, 
	this
	};
    DataHandle<TrajectoriesContainer> m_inputTrajectories {
      "inputTrajectories", 
	Gaudi::DataHandle::Reader, 
	this
	};

    // Output collection;
    std::unique_ptr<DataHandle<eic::TofParticleIDCollection>> m_outputTofPID;

    /// Pointer to the geometry and PDG service;
    SmartIF<IGeoSvc> m_geoSvc;
    SmartIF<IParticleSvc> m_pidSvc;

    // A .root file with a CherenkovDetectorCollection entry, to run in a back door mode;
    Gaudi::Property<std::string>            m_Detector            {this, "Detector",             "BarrelTOF"};
    // Be careful with the units here!;
    Gaudi::Property<double>                 m_Threshold           {this, "Threshold",   0 * Gaudi::Units::eV};
    Gaudi::Property<double>                 m_Resolution          {this, "Resolution",  0 * Gaudi::Units::ps};

  private:
    Rndm::Numbers m_rngGauss, m_rngUni;

    Acts::GeometryContext m_geoContext;
    Acts::MagneticFieldContext m_fieldContext;

    BoundTrackParamPtrResult propagateTrack(const Acts::BoundTrackParameters& params, 
					    const SurfacePtr& targetSurf);
  };

  DECLARE_COMPONENT(TOFPID)

} // namespace Jug::PID

#endif
