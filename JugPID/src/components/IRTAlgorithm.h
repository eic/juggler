
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/IParticleSvc.h"
#include "JugBase/UniqueID.h"

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/PhotoMultiplierHitCollection.h"

#include "eicd/TrajectoryCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/CherenkovParticleIDCollection.h"

//#define _USE_RECONSTRUCTED_TRACKS_
//#define _USE_TRAJECTORIES_

#ifndef _IRT_ALGORITHM_
#define _IRT_ALGORITHM_

#include "IRTAlgorithmServices.h"

class CherenkovRadiator;
class CherenkovDetector;
class CherenkovDetectorCollection;

//typedef SimpleProperty   < std::vector< std::pair<std::string, double> > >      QStringArrayProperty;
//typedef SimplePropertyRef< std::vector< std::pair<std::string, double> > >      QStringArrayPropertyRef;

namespace Jug::PID {
  /** RICH IRT
   *
   * \ingroup reco
   */
  class IRTAlgorithm : public IRTAlgorithmServices, public GaudiAlgorithm, AlgorithmIDMixin<> {

  public:
    // constructor
    IRTAlgorithm(const std::string& name, ISvcLocator* svcLoc);

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
#ifdef _USE_TRAJECTORIES_
    DataHandle<eic::TrajectoryCollection> m_inputTrajectories {
      "inputTrajectories", 
	Gaudi::DataHandle::Reader, 
	this
	};
#endif

    // Output collection;
    std::unique_ptr<DataHandle<eic::CherenkovParticleIDCollection>> m_outputCherenkovPID;

    /// Pointer to the geometry and PDG service;
    SmartIF<IGeoSvc> m_geoSvc;
    SmartIF<IParticleSvc> m_pidSvc;

    // A .root file with a CherenkovDetectorCollection entry, to run in a back door mode;
    Gaudi::Property<std::string>                            m_ConfigFile          {this, "ConfigFile", ""};
    Gaudi::Property<std::string>                            m_Detector            {this, "Detector", ""};
    // The famous ~0.85 for the S13660-3050AE-08 SiPM and its QE table expected;
    Gaudi::Property<double>                                 m_GeometricEfficiency {this, "GeometricEfficiency", 1.0};
    // Array of {e, qe} points and the desired number of bins in a fast lookup table;
    Gaudi::Property<std::vector<std::pair<double, double>>> m_QE_input_data       {this, "QEcurve"};
    Gaudi::Property<unsigned>                               m_QEbins              {this, "QEbins",      0};
    Gaudi::Property<std::vector<std::string>>               m_RadiatorConfig      {this, "Radiators"};

  private:
    // FIXME: will need several detectors;
    CherenkovDetectorCollection *m_IrtGeo;
    CherenkovDetector           *m_IrtDet;

    std::vector<CherenkovRadiator*> m_SelectedRadiators;

    Rndm::Numbers m_rngUni;
  };

  DECLARE_COMPONENT(IRTAlgorithm)

} // namespace Jug::PID

#endif
