
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/IParticleSvc.h"
#include "JugBase/UniqueID.h"

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/PhotoMultiplierHitCollection.h"

#include "eicd/ReconstructedParticleCollection.h"
//#include "eicd/CherenkovParticleIDCollection.h"

//#define _USE_RECONSTRUCTED_TRACKS_

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
    std::unique_ptr<DataHandle<dd4pod::PhotoMultiplierHitCollection>> m_inputHitCollection;
    DataHandle<dd4pod::Geant4ParticleCollection> m_inputMCParticles {
      "inputMCParticles", 
	Gaudi::DataHandle::Reader, 
	this
	};
#ifdef _USE_RECONSTRUCTED_TRACKS_
    DataHandle<eic::ReconstructedParticleCollection> m_inputRecoParticles {
      "inputRecoParticles", 
	Gaudi::DataHandle::Reader, 
	this
	};
#endif

    // Output collection;
    //+std::unique_ptr<DataHandle<eic::CherenkovParticleIDCollection>> m_outputCherenkovPID;

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

  };

  DECLARE_COMPONENT(TOFPID)

} // namespace Jug::PID

#endif
