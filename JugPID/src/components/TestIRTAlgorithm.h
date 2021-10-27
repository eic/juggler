
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/IParticleSvc.h"
#include "JugBase/UniqueID.h"

using namespace Gaudi::Units;

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/PhotoMultiplierHitCollection.h"

#include "eicd/TrajectoryCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/CherenkovParticleIDCollection.h"

#define _USE_RECONSTRUCTED_TRACKS_

#ifndef _TEST_IRT_ALGORITHM_
#define _TEST_IRT_ALGORITHM_

#include "TestIRTAlgorithmServices.h"

class CherenkovDetector;
class CherenkovDetectorCollection;

namespace Jug::Reco {
  /** RICH IRT
   *
   * \ingroup reco
   */
  class TestIRTAlgorithm : public TestIRTAlgorithmServices, public GaudiAlgorithm, AlgorithmIDMixin<> {

  public:
    // constructor
    TestIRTAlgorithm(const std::string& name, ISvcLocator* svcLoc);

    // INITIALIZE .....................................................
    StatusCode initialize( void ) override;

    // EXECUTE .....................................................
    StatusCode execute( void ) override;

    // FINALIZE .....................................................
    StatusCode finalize( void ) override;

    // FIXME: the collection names make no sense; synchronize with the upstream plugins;
    // FIXME: names either here or in the source code must be duplicates;
    DataHandle<dd4pod::PhotoMultiplierHitCollection> m_inputHitCollection {
      "inputHitCollection",
	Gaudi::DataHandle::Reader,
	this
	};
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
    DataHandle<eic::TrajectoryCollection> m_inputTrajectories {
      "inputTrajectories", 
	Gaudi::DataHandle::Reader, 
	this
	};

    DataHandle<eic::CherenkovParticleIDCollection> m_outputCherenkovPID {
      "CherenkovPID", 
	Gaudi::DataHandle::Writer, 
	this
	};

    /// Pointer to the geometry and PDG service;
    SmartIF<IGeoSvc> m_geoSvc;
    SmartIF<IParticleSvc> m_pidSvc;

    // A .root file with a CherenkovDetectorCollection entry, to run in a back door mode;
    Gaudi::Property<std::string>                            m_ConfigFile          {this, "ConfigFile", ""};
    // The famous ~0.85 for the S13660-3050AE-08 SiPM and its QE table expected;
    Gaudi::Property<double>                                 m_GeometricEfficiency {this, "GeometricEfficiency", 1.0};
    // Array of {e, qe} points and the desired number of bins in a fast lookup table;
    Gaudi::Property<std::vector<std::pair<double, double>>> m_QE_input_data       {this, "QEcurve"};
    Gaudi::Property<unsigned>                               m_QEbins              {this, "QEbins",      0};

  private:
    // FIXME: will need several detectors;
    CherenkovDetectorCollection *m_IrtGeo;
    CherenkovDetector           *m_IrtDet;

    Rndm::Numbers m_rngUni;//, m_rngNorm;
  };

  DECLARE_COMPONENT(TestIRTAlgorithm)

} // namespace Jug::Reco

#endif
