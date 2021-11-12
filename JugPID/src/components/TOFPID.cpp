
#include <map>

#include <TVector3.h>

#include "DDRec/CellIDPositionConverter.h"

#include "TOFPID.h"

// FIXME: make sure that 'mm' in initialize() are what they are expected!;
using namespace Gaudi::Units;

// -------------------------------------------------------------------------------------

Jug::PID::TOFPID::TOFPID(const std::string& name, ISvcLocator* svcLoc) 
  : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info())
{
  declareProperty("inputMCParticles",         m_inputMCParticles,              "");
  declareProperty("inputRecoParticles",       m_inputRecoParticles,            "");
  declareProperty("inputTrajectories",        m_inputTrajectories,             "");
} // Jug::PID::TOFPID::TOFPID()

// -------------------------------------------------------------------------------------

BoundTrackParamPtrResult Jug::PID::TOFPID::propagateTrack(const Acts::BoundTrackParameters& params, 
							  const SurfacePtr& targetSurf) 
{ 
#if _OFF_    
  std::cout << "Propagating final track fit with momentum: " 
            << params.momentum() << " and position " 
            << params.position(m_geoContext)
            << std::endl
            << "track fit phi/eta "
            << atan2(params.momentum()(1), 
		     params.momentum()(0)) 
            << " and " 
            << atanh(params.momentum()(2) 
		     / params.momentum().norm())
            << std::endl;
#endif

  //+@@@std::shared_ptr<const Acts::TrackingGeometry>
  //+@@@      trackingGeometry = m_geoSvc->trackingGeometry();
  std::shared_ptr<const Acts::MagneticFieldProvider>
    magneticField = m_geoSvc->getFieldProvider();
  
  using Stepper            = Acts::EigenStepper<>;
  using Propagator         = Acts::Propagator<Stepper>;
  
  Stepper stepper(magneticField);
  Propagator propagator(stepper);
  
  // Acts::Logging::Level logLevel = Acts::Logging::FATAL;
  Acts::Logging::Level logLevel = Acts::Logging::DEBUG;
  
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("ProjectTrack Logger", logLevel));
  
  Acts::PropagatorOptions<> options(m_geoContext,
				    m_fieldContext,
				    Acts::LoggerWrapper{logger()});
  
  auto result = propagator.propagate(params, *targetSurf, options);
  
  if(result.ok()) return std::move((*result).endParameters);
  
  return result.error();
} // Jug::PID::TOFPID::propagateTrack()

// -------------------------------------------------------------------------------------

StatusCode Jug::PID::TOFPID::initialize( void ) 
{
  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  } //if

  // Geometry service, to retrieve materials, locations, etc; see e.g. IRTAlgorithm.cpp;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  // Particle service to retries charges, masses, etc; 
  m_pidSvc = service("ParticleSvc");
  if (!m_pidSvc) {
    error() << "Unable to locate Particle Service. "
            << "Make sure you have ParticleSvc in the configuration."
            << endmsg;
    return StatusCode::FAILURE;
  } //if

  // Need a random number generator too, perhaps;
  {
    auto randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);

    // FIXME: check Rndm::Gauss() parameter order;
    auto gauss = m_rngGauss.initialize(randSvc, Rndm::Gauss(0., 1.));
    auto flat  = m_rngUni.initialize  (randSvc, Rndm::Flat (0., 1.));
    if (!gauss.isSuccess() || !flat.isSuccess()) {
      error() << "Cannot initialize random generator!" << endmsg;
      return StatusCode::FAILURE;
    } //if
  }

  // Detector name ("BarrelTOF" presumably);
  std::string dname = m_Detector.value();

  // Input hit collection; detector name tagged;
  m_inputHitCollection =
    std::make_unique<DataHandle<dd4pod::TrackerHitCollection>>((dname + "Hits").c_str(), 
							       Gaudi::DataHandle::Reader, this);
  // Output PID info collection;
  m_outputTofPID =
    std::make_unique<DataHandle<eic::TofParticleIDCollection>>((dname + "PID").c_str(), 
							       Gaudi::DataHandle::Writer, this);

  return StatusCode::SUCCESS;
} // Jug::PID::TOFPID::initialize()

// -------------------------------------------------------------------------------------

StatusCode Jug::PID::TOFPID::execute( void )
{
  // Input collection(s);
  const auto &hits         = *m_inputHitCollection->get();
  const auto &mctracks     = *m_inputMCParticles.get();
  const auto &rctracks     = *m_inputRecoParticles.get();
  const TrajectoriesContainer* trajectories = m_inputTrajectories.get();
  // Output collection(s);
  auto &cpid               = *m_outputTofPID->createAndPut();
    
  // Loop through the reconstructed tracks in this code; obviously only for them
  // we may have ACTS trajectory parameterization; 
  printf("@@@ tracks: %ld\n", rctracks.size());
  printf("@@@ hits: %ld\n",   hits.size());
  for(const auto &rctrack: rctracks) {
    unsigned id = rctrack.mcID().value;
    // FIXME: do the sanity checks in a more meaningful way later;
    if (id >= mctracks.size()) continue;
    
    // Link to the MC track;
    const auto &mctrack = mctracks[id];

    // Now just follow the logic of TrackParamTruthInit.cpp; 

    // genStatus = 1 means thrown G4Primary, but dd4gun uses genStatus == 0
    if (mctrack.genStatus() > 1 ) {
      if (msgLevel(MSG::DEBUG)) {
	debug() << "ignoring particle with genStatus = " << mctrack.genStatus() << endmsg;
      }
      continue;
    } //if

    // This is here just to show how to use m_pidSvc;
    const double charge = m_pidSvc->particle(mctrack.pdgID()).charge;
    if (abs(charge) < std::numeric_limits<double>::epsilon()) {
      if (msgLevel(MSG::DEBUG)) {
	debug() << "ignoring neutral particle" << endmsg;
      }
      continue;
    } //if

      // FIXME: consider only primaries for the time being;
    if (mctrack.g4Parent()) continue;

    //printf("@@@ %d\n", rctrack.ID().value);

    // Get track parameterization at a given radius; Zhenyu, I'm not sure indexing 
    // correspondence here is a 1:1; perhaps a remapping will be required; start with the 
    // single-track events;
    const auto& traj = (*trajectories)[rctrack.ID().value];
    const auto& mj        = traj.multiTrajectory();
    const auto& trackTips = traj.tips();
    if (trackTips.empty()) {
      printf("@@@ trackTips.empty(): should not happen?\n");
      continue;
    } //if

    //printf("@@@ %ld\n", trackTips.size());
    auto& trackTip = trackTips.front();

    if (traj.hasTrackParameters(trackTip)) {
      const auto& boundParam = traj.trackParameters(trackTip);
      //const auto& parameter  = boundParam.parameters();
      //debug() << "pos 0 = " << parameter[Acts::eFreePos0] << endmsg;
      //debug() << "pos 1 = " << parameter[Acts::eFreePos1] << endmsg;
      //debug() << "pos 2 = " << parameter[Acts::eFreePos2] << endmsg;
	    
      // make a reference cylinder to mimic BTOF; one can use m_geoSvc to extract this number;
      const auto btofRadius = 500.0;
      const auto btofEta = 1.1;
      const auto btofTheta = 2. * atan(exp(-btofEta));
      const auto halfZ = btofRadius / tan(btofTheta);  
      std::shared_ptr<Acts::DiscSurface> btofSurf = 
	Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), btofRadius, halfZ);
      
      auto result = propagateTrack(boundParam, btofSurf);
      if(result.ok()) {
	auto trackStateParams = std::move(**result);
	auto projectionPos = trackStateParams.position(m_geoContext);
	auto projectionMom = trackStateParams.momentum();//m_geoContext);
	
	// And one also needs to fish out the trajectory path length to this point;
	TVector3 xx (projectionPos(0), projectionPos(1), projectionPos(2));
	TVector3 p0 (projectionMom(0), projectionMom(1), projectionMom(2));
	TVector3 nn = p0.Unit();
      } else {
	printf("@@@ result.ok() = 0: take care about this possibility\n");
	continue;
      } //if
    } //if

    // Now find a matching (Barrel) TOF hit, perform TOF PID calculations, and fill 
    // out the output structures;
    for(const auto &hit: hits) {
      // FIXME: yes, use MC truth here; FIXME: should create a lookup std::map;
      if (hit.g4ID() != mctrack.ID()) continue;
      
      {
	const auto &x = hit.position();
	auto hid = hit.cellID();
      }

      {
	// Create an output PID entry;
	auto cbuffer = cpid.create();

	// Some random number; can be a confidence level for the "winning hypothesis", etc;
	cbuffer.reserved(111.1);

	// Fill out just one single (pion) hypothesis as an example; for more than one 
	// there should just be a loop with more than one cbuffer.addoptions(hypothesis) call;
	eic::TofPdgHypothesis hypothesis;
	
	// Let's say '1' stands for the barrel TOF;
	hypothesis.detector = 1;
	// pion;
	hypothesis.pdg      = 211;
	// Something like a chi^ CDF, which come after the track-to-hit matching in time
	// for a given mass hypothesis;
	hypothesis.chi2cdf  = 0.55;
	// Some random number; 
	hypothesis.reserved = 222.2;
	
	cbuffer.addoptions(hypothesis);
	
	// FIXME: well, and what does go here instead of 0?;
	cbuffer.ID({0, algorithmID()});
	
	// Reference to the reconstructed track;
	cbuffer.recID(rctrack.ID());
      }  
    } //for hit
  } //for rctrack

  return StatusCode::SUCCESS;
} // Jug::PID::TOFPID::execute()

// -------------------------------------------------------------------------------------

StatusCode Jug::PID::TOFPID::finalize( void ) 
{
  info() << "TOFPID: Finalizing..." << endmsg;

  return Algorithm::finalize(); // must be executed last
} // Jug::PID::TOFPID::finalize()

// -------------------------------------------------------------------------------------


