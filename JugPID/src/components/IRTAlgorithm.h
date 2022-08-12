
//#define _USE_RECONSTRUCTED_TRACKS_
//#define _USE_STORED_TRAJECTORIES_
//#define _USE_ON_THE_FLY_TRAJECTORIES_

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/IParticleSvc.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "eicd/TrajectoryCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/CherenkovParticleIDCollection.h"

#ifdef _USE_ON_THE_FLY_TRAJECTORIES_
#include "JugBase/BField/DD4hepBField.h"
#include "JugTrack/Trajectories.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"

using BoundTrackParamPtr = std::unique_ptr<const Acts::BoundTrackParameters>;
using BoundTrackParamPtrResult = Acts::Result<BoundTrackParamPtr>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
#endif

#ifndef _IRT_ALGORITHM_
#define _IRT_ALGORITHM_

#include "IRTAlgorithmServices.h"

class CherenkovRadiator;
class CherenkovDetector;
class CherenkovDetectorCollection;

namespace Jug::PID {
  /** RICH IRT
   *
   * \ingroup reco
   */
  class IRTAlgorithm : public IRTAlgorithmServices, public GaudiAlgorithm {

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
    std::unique_ptr<DataHandle<edm4hep::SimTrackerHitCollection>> m_inputHitCollection;
    DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles {
      "inputMCParticles", 
	Gaudi::DataHandle::Reader, 
	this
	};
#ifdef _USE_RECONSTRUCTED_TRACKS_
    DataHandle<eicd::ReconstructedParticleCollection> m_inputRecoParticles {
      "inputRecoParticles", 
	Gaudi::DataHandle::Reader, 
	this
	};
#endif
#ifdef _USE_STORED_TRAJECTORIES_
    DataHandle<eicd::TrajectoryCollection> m_inputTrajectories {
      "inputTrajectories", 
	Gaudi::DataHandle::Reader, 
	this
	};
#endif
#ifdef _USE_ON_THE_FLY_TRAJECTORIES_
    DataHandle<TrajectoriesContainer> m_inputTrajectories {
      "inputTrajectories", 
	Gaudi::DataHandle::Reader, 
	this
	};
#endif

    // Output collection;
    std::unique_ptr<DataHandle<eicd::CherenkovParticleIDCollection>> m_outputCherenkovPID;

    /// Pointer to the geometry and PDG service;
    SmartIF<IGeoSvc> m_geoSvc;
    SmartIF<IParticleSvc> m_pidSvc;

    // A .root file with a CherenkovDetectorCollection entry, to run in a back door mode;
    Gaudi::Property<std::string>                            m_ConfigFile          {this, "ConfigFile", ""};
    Gaudi::Property<std::string>                            m_Detector            {this, "Detector", ""};
    // The famous ~0.85 for the S13660-3050AE-08 SiPM and its QE table expected; well, as long 
    // as the fill factor is actually mimiced in E(D)Rich_geo.cpp, geometric efficiency should be 1.0 here;
    //Gaudi::Property<double>                                 m_GeometricEfficiency {this, "GeometricEfficiency", 1.0};
    // Another famous factor (life is always below expectations) which people like to see ~0.7;
    Gaudi::Property<double>                                 m_SafetyFactor        {this, "SafetyFactor", 1.0};
    Gaudi::Property<double>                                 m_SensorPixelSize     {this, "SensorPixelSize",  3.0};
    Gaudi::Property<double>                                 m_SensorPixelPitch    {this, "SensorPixelPitch", 3.2};
    // Array of {e, qe} points and the desired number of bins in a fast lookup table;
    Gaudi::Property<std::vector<std::pair<double, double>>> m_QE_input_data       {this, "QEcurve"};
    Gaudi::Property<unsigned>                               m_QEbins              {this, "QEbins",      0};
    Gaudi::Property<std::vector<std::string>>               m_RadiatorConfig      {this, "Radiators"};

  private:
    // FIXME: 'const' does not seem to work (?), since this structure is populated via the front 
    // door mechanism, piece by piece, instead of being assigned once;
    /*const*/ CherenkovDetectorCollection *m_IrtGeo;
    /*const*/ CherenkovDetector           *m_IrtDet;

    // This is needed for dd4hep cell index decoding;
    uint64_t m_ReadoutCellMask;

    std::vector<CherenkovRadiator*> m_SelectedRadiators;

    Rndm::Numbers m_rngGauss, m_rngUni;

#ifdef _USE_ON_THE_FLY_TRAJECTORIES_
    Acts::GeometryContext m_geoContext;
    Acts::MagneticFieldContext m_fieldContext;

    BoundTrackParamPtrResult propagateTrack(const Acts::BoundTrackParameters& params, 
					    const SurfacePtr& targetSurf);
#endif

  };

  DECLARE_COMPONENT(IRTAlgorithm)

} // namespace Jug::PID

#endif
