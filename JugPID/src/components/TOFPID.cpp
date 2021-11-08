
#include <map>

#include "DDRec/CellIDPositionConverter.h"

#include "TOFPID.h"

// FIXME: make sure that 'mm' in initialize() are what they are expected!;
using namespace Gaudi::Units;

// -------------------------------------------------------------------------------------

Jug::PID::TOFPID::TOFPID(const std::string& name, ISvcLocator* svcLoc) 
  : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info())
{
  declareProperty("inputMCParticles",                 m_inputMCParticles,              "");
#ifdef _USE_RECONSTRUCTED_TRACKS_
  declareProperty("inputRecoParticles",               m_inputRecoParticles,            "");
#endif
} // Jug::PID::TOFPID::TOFPID()

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
    std::make_unique<DataHandle<dd4pod::PhotoMultiplierHitCollection>>((dname + "Hits").c_str(), 
									 Gaudi::DataHandle::Reader, this);
  // Output PID info collection;
  //m_outputCherenkovPID =
  //std::make_unique<DataHandle<eic::CherenkovParticleIDCollection>>((dname + "PID").c_str(), 
  //								     Gaudi::DataHandle::Writer, this);

  return StatusCode::SUCCESS;
} // Jug::PID::TOFPID::initialize()

// -------------------------------------------------------------------------------------

StatusCode Jug::PID::TOFPID::execute( void )
{
  // Input collection(s);
  const auto &hits         = *m_inputHitCollection->get();
  const auto &mctracks     = *m_inputMCParticles.get();
#ifdef _USE_RECONSTRUCTED_TRACKS_
  const auto &rctracks     = *m_inputRecoParticles.get();
#endif
  // Output collection(s);
  //+auto &cpid               = *m_outputCherenkovPID->createAndPut();
    
  // Loop through either MC or reconstructed tracks;
#ifdef _USE_RECONSTRUCTED_TRACKS_
  for(const auto &rctrack: rctracks) {
    unsigned id = rctrack.mcID().value;
    // FIXME: do the sanity checks in a more meaningful way later;
    if (id >= mctracks.size()) continue;
    
    const auto &mctrack = mctracks[id];
#else
  for(const auto &mctrack: mctracks) {
#endif

    // At this point have a reference to a 'mctrack', either this or that way;
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

    for(const auto &hit: hits) {
      // FIXME: yes, use MC truth here; not really needed I guess; 
      if (hit.g4ID() != mctrack.ID()) continue;
      
      // Simulate QE & geometric sensor efficiency; FIXME: hit.energy() is numerically 
      // in GeV units, but Gaudi::Units::GeV = 1000; prefer to convert photon energies 
      // to [eV] in all places by hand;
      //if (!QE_pass(1E9*hit.energy(), m_rngUni()*m_GeometricEfficiency.value())) 
      //continue;
      
      {
	//const auto &x = hit.position();
	//auto id = hit.cellID();
      }
    } //for hit

#if _TODAY_
    {
      auto cbuffer = cpid.create();
      eic::CherenkovPdgHypothesis hypothesis;
	  	  
      hypothesis.npe    = hypo->GetNpe   (radiator);
      hypothesis.weight = hypo->GetWeight(radiator);
	  
      cbuffer.addoptions(hypothesis);
	
      // FIXME: well, and what does go here instead of 0?;
      cbuffer.ID({0, algorithmID()});
	
	// Reference to either MC track or a reconstructed track;
#ifdef _USE_RECONSTRUCTED_TRACKS_
      cbuffer.recID(rctrack.ID());
#else
      cbuffer.recID(mctrack.ID());
#endif
    }  
#endif
  } //for rctrack (mctrack)

  return StatusCode::SUCCESS;
} // Jug::PID::TOFPID::execute()

// -------------------------------------------------------------------------------------

StatusCode Jug::PID::TOFPID::finalize( void ) 
{
  info() << "TOFPID: Finalizing..." << endmsg;

  return Algorithm::finalize(); // must be executed last
} // Jug::PID::TOFPID::finalize()

// -------------------------------------------------------------------------------------


