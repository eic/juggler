
#include <map>

#include "TFile.h"

#include "DDRec/CellIDPositionConverter.h"

#include "IRTAlgorithm.h"

#include "IRT/CherenkovRadiator.h"
#include "IRT/CherenkovEvent.h"
#include "IRT/CherenkovDetectorCollection.h"

// FIXME: make sure that 'mm' in initialize() are what they are expected!;
using namespace Gaudi::Units;

// -------------------------------------------------------------------------------------

Jug::PID::IRTAlgorithm::IRTAlgorithm(const std::string& name, ISvcLocator* svcLoc) 
  : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()), m_IrtGeo(0), m_IrtDet(0), 
    m_ReadoutCellMask(0x0)
{
  declareProperty("inputMCParticles",                 m_inputMCParticles,              "");
#ifdef _USE_RECONSTRUCTED_TRACKS_
  declareProperty("inputRecoParticles",               m_inputRecoParticles,            "");
#endif
#ifdef _USE_STORED_TRAJECTORIES_
  declareProperty("inputTrajectories",                m_inputTrajectories,             "");
#endif
} // Jug::PID::IRTAlgorithm::IRTAlgorithm()

// -------------------------------------------------------------------------------------

//
// FIXME: several things need to be cached here once they are made working;
//

#ifdef _USE_ON_THE_FLY_TRAJECTORIES_
BoundTrackParamPtrResult Jug::PID::IRTAlgorithm::propagateTrack(const Acts::BoundTrackParameters& params, 
								const SurfacePtr& targetSurf) 
{ 
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
  
  auto result = propagator.propagate(params, *targetSurf, 
				     options);
  
  if(result.ok()) return std::move((*result).endParameters);
  
  return result.error();
} // Jug::PID::IRTAlgorithm::propagateTrack()
#endif
// -------------------------------------------------------------------------------------

StatusCode Jug::PID::IRTAlgorithm::initialize( void ) 
{
  // Either get rid of some of this branching, after the most reasonable way is 
  // chosen, or just make the #define's configurable;
#if defined (_USE_STORED_TRAJECTORIES_) && defined(_USE_ON_THE_FLY_TRAJECTORIES_)
#error _USE_STORED_TRAJECTORIES_ and _USE_ON_THE_FLY_TRAJECTORIES_ are mutually exclusive.
#endif
#if defined(_USE_ON_THE_FLY_TRAJECTORIES_) && !defined(_USE_RECONSTRUCTED_TRACKS_)
#error _USE_ON_THE_FLY_TRAJECTORIES_ requires _USE_RECONSTRUCTED_TRACKS_.
#endif

  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  } //if

  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  
  m_pidSvc = service("ParticleSvc");
  if (!m_pidSvc) {
    error() << "Unable to locate Particle Service. "
            << "Make sure you have ParticleSvc in the configuration."
            << endmsg;
    return StatusCode::FAILURE;
  } //if

  {
    std::string config = m_ConfigFile.value();
    std::string dname = m_Detector.value();

    // Well, a back door: if a config file name was given, import from file;
    if (config.size()) {
      auto fcfg  = new TFile(config.c_str());
      if (!fcfg) {
        error() << "Failed to open IRT config .root file." << endmsg;
        return StatusCode::FAILURE;
      } //if

      m_IrtGeo = dynamic_cast<CherenkovDetectorCollection*>(fcfg->Get("CherenkovDetectorCollection"));
      if (!m_IrtGeo) {
        error() << "Failed to import IRT geometry from the config .root file." << endmsg;
        return StatusCode::FAILURE;
      } //if

      // May want to fish the detector name out of the geometry file, if it is unique;
      if (dname.empty()) {
        auto &detectors = m_IrtGeo->GetDetectors();
        if (detectors.size() != 1) {
          error() << "More than one detector in the provided IRT geometry config .root file." << endmsg;
          return StatusCode::FAILURE;
        } //if

        dname = (*detectors.begin()).first.Data();
        m_IrtDet = (*detectors.begin()).second;
      } else {
        // Detector pointer;
        m_IrtDet = m_IrtGeo->GetDetector(dname.c_str());
        if (!m_IrtDet) {
          error() << "Failed to import IRT geometry from the config .root file." << endmsg;
          return StatusCode::FAILURE;
        } //if
      } //if

      m_ReadoutCellMask = m_IrtDet->GetReadoutCellMask();
    } //if

    if (dname.empty()) {
      error() << "No RICH detector name provided." << endmsg;
      return StatusCode::FAILURE;
    } //if

    {
      // Extract refractive index information; build internal lookup table;
      auto dd4Det = m_geoSvc->detector();
      
      for(auto rptr: m_IrtDet->Radiators()) {
	const auto radiator = rptr.second;
	// FIXME: sanity checks;
	auto rindex = dd4Det->material(radiator->GetAlternativeMaterialName()).property("RINDEX"); 
	std::vector<std::pair<double, double>> buffer;
	for(unsigned iq=0; iq<rindex->GetRows(); iq++)
	  buffer.push_back(std::make_pair(1E9*rindex->Get(iq,0), rindex->Get(iq,1)));

	// FIXME: may want to use different binning here?;
	radiator->m_ri_lookup_table = ApplyFineBinning(buffer, m_QEbins.value());
      } //for rptr
    } 
	  
    // Input hit collection;
    m_inputHitCollection =
      std::make_unique<DataHandle<dd4pod::PhotoMultiplierHitCollection>>((dname + "Hits").c_str(), 
									 Gaudi::DataHandle::Reader, this);
    // Output PID info collection;
    m_outputCherenkovPID =
      std::make_unique<DataHandle<eic::CherenkovParticleIDCollection>>((dname + "PID").c_str(), 
								       Gaudi::DataHandle::Writer, this);
  }

    // FIXME: how about Gaudi::Property<std::vector<std::tuple>>?;
  {
    const auto &radiators = m_RadiatorConfig.value();
    
    // radiator specifications, to be parsed from options file
    std::string name = "";
    std::string mode = "";
    unsigned zbins = 0;
    float qmrad = -1.;
    float ri = -1.;
    
    for(auto rptr: radiators) {
      // parse radiator specification string; FIXME: expected format could be more consistent, e.g. each token is `this=that`
      //                                      FIXME: check for incorrect formatting
      /* expected tokens (delimited by space, order does not matter):
       * - "%s"          parses to  `name`
       * - "zbins=%u"    parses to  `zbins`
       * - "smearing=%s" parses to  `mode`
       * - "%fmrad"      parses to  `mrad`
       */
      std::stringstream rptrStream(rptr);
      std::string rptrTok;
      std::regex rmEq("^.*=");
      while(rptrStream >> rptrTok) {
	if(rptrTok.find("zbins")!=std::string::npos) 
	  zbins = std::stoul(std::regex_replace(rptrTok,rmEq,""));
	else if(rptrTok.find("rindex")!=std::string::npos) 
	  ri = std::stof(std::regex_replace(rptrTok,rmEq,""));
	else if(rptrTok.find("smearing")!=std::string::npos) 
	  mode = std::regex_replace(rptrTok,rmEq,"");
	else if(rptrTok.find("mrad")!=std::string::npos) 
	  qmrad = std::stof(std::regex_replace(rptrTok,std::regex("mrad$"),""));
	else name = rptrTok; // otherwise assume it's the name
      }
      //printf("@@@ %s %s %d %7.1f\n", name.c_str(), mode.c_str(), zbins, qmrad);
      
      if (!zbins || qmrad<0.0 || ri<0.0 || name.compare("")==0 || mode.compare("")==0) {
	error() << "Parsed radiator parameters do not pass a sanity check." << endmsg;
	return StatusCode::FAILURE;
      } //if
      
	// FIXME: sanity check would not hurt;
      auto radiator = m_IrtDet->GetRadiator(name.c_str());
      if (radiator) {
	m_SelectedRadiators.push_back(radiator);
	
	// FIXME: this looks awkward;
	radiator->m_AverageRefractiveIndex = ri;//radiator->n();
	radiator->SetReferenceRefractiveIndex(ri);

	if (qmrad) {
	  // Well, want them readable in the config file, but pass in [rad] to the IRT algorithm;
	  qmrad /= 1000.;
	  
	  if (strcmp(mode.c_str(), "uniform") && strcmp(mode.c_str(), "gaussian")) {
	    error() << "Unknown smearing mode." << endmsg;
	    return StatusCode::FAILURE;
	  } //if
	  
	  (!strcmp(mode.c_str(), "uniform") ? radiator->SetUniformSmearing(qmrad) : radiator->SetGaussianSmearing(qmrad));
	} //if
	
	radiator->SetTrajectoryBinCount(zbins);
      } else {
	error() << "Unknown radiator name." << endmsg;
	return StatusCode::FAILURE;
      } //if
    } //for radiator
    
    if (!m_SelectedRadiators.size()) {
      error() << "At least one radiator name should be provided." << endmsg;
      return StatusCode::FAILURE;
    } //if
  }

  // Need a random number generator for the QE stuff;
  {
    auto randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
    auto gauss = m_rngGauss.initialize(randSvc, Rndm::Gauss(0., 1.));     
    auto flat  = m_rngUni.initialize  (randSvc, Rndm::Flat (0., 1.));     
    if (!gauss.isSuccess() || !flat.isSuccess()) {       
      error() << "Cannot initialize random generator!" << endmsg;       
      return StatusCode::FAILURE;     
    } //if
  }

  // Initialize the QE lookup table;
  m_QE_lookup_table = ApplyFineBinning(m_QE_input_data.value(), m_QEbins.value());

  return StatusCode::SUCCESS;
} // Jug::PID::IRTAlgorithm::initialize()

// -------------------------------------------------------------------------------------

StatusCode Jug::PID::IRTAlgorithm::execute( void )
{
  // Input collection(s);
  const auto &hits         = *m_inputHitCollection->get();
  const auto &mctracks     = *m_inputMCParticles.get();
#ifdef _USE_RECONSTRUCTED_TRACKS_
  const auto &rctracks     = *m_inputRecoParticles.get();
#endif
#ifdef _USE_STORED_TRAJECTORIES_
  const auto &trajectories = *m_inputTrajectories.get();
  // First populate the trajectory-to-reconstructed (or -to-simulated) mapping table;
  std::map<eic::Index, const eic::ConstTrajectory*> rc2trajectory;
  for(const auto &trajectory: trajectories)
    rc2trajectory[trajectory.trackID()] = &trajectory;
#endif
#ifdef _USE_ON_THE_FLY_TRAJECTORIES_
  const TrajectoriesContainer* trajectories = m_inputTrajectories.get();
#endif
  // Output collection(s);
  auto &cpid               = *m_outputCherenkovPID->createAndPut();
    
  // An interface variable; FIXME: check memory cleanup;
  auto event = new CherenkovEvent();

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

    const double charge = m_pidSvc->particle(mctrack.pdgID()).charge;
    if (abs(charge) < std::numeric_limits<double>::epsilon()) {
      if (msgLevel(MSG::DEBUG)) {
	debug() << "ignoring neutral particle" << endmsg;
      }
      continue;
    } //if

      // FIXME: consider only primaries for the time being;
    if (mctrack.g4Parent()) continue;
    //printf("@@@ %d %d %d\n", mctrack.ID(), mctrack.pdgID(), mctrack.g4Parent());

    auto particle = new ChargedParticle(mctrack.pdgID());
    event->AddChargedParticle(particle);

    // Start radiator history for all known radiators;
    for(auto rptr: m_IrtDet->Radiators()) 
      particle->StartRadiatorHistory(std::make_pair(rptr.second, new RadiatorHistory()));
      
    // Distribute photons over radiators where they were presumably produced, to 
    // create a structure which would mimic a standalone G4 stepper behavior;
    bool useful_photons_found = false;
    for(const auto &hit: hits) {
      // FIXME: range checks;
      const auto &phtrack = mctracks[hit.truth().trackID];

      // FIXME: yes, use MC truth here; not really needed I guess;
      // FIXME: range checks;
      if (phtrack.parents()[0] != mctrack.ID()) continue;
      
      // Vertex where photon was created;
      const auto &vs = phtrack.vs();
      TVector3 vtx(vs.x, vs.y, vs.z);
      
      // FIXME: this is in general correct, but needs refinement;
      TVector3 ip(0,0,0);
      // FIXME: unify with the code below;
      auto radiator = m_IrtDet->GuessRadiator(vtx, (vtx - ip).Unit());
      // FIXME: do it better later;
      if (!radiator) continue;

      // Simulate QE & geometric sensor efficiency; FIXME: hit.energy() is numerically 
      // in GeV units, but Gaudi::Units::GeV = 1000; prefer to convert photon energies 
      // to [eV] in all places by hand;
      double eVenergy = 1E9*hit.energy();
      if (!QE_pass(eVenergy, m_rngUni()) || 
	  m_rngUni() > m_GeometricEfficiency.value()*m_SafetyFactor.value()) {

	// Prefer to create the stepping history in a clean way, namely use only those 
	// photons which were rejected; FIXME: optionally include all; anyway, now 
	// figure out (1) to which detector sector does this 3D point belongs, (2) 
	// which radiator this could have been;

	// Add the point to the radiator history buffer; FIXME: all this is not 
	// needed if trajectory parameterization is used;
	particle->FindRadiatorHistory(radiator)->AddStepBufferPoint(phtrack.time(), vtx);

	continue;
      } //if

      useful_photons_found = true;

      // Photon accepted; add it to the internal event structure;
      auto photon = new OpticalPhoton();
      
      // FIXME: '0' - for now assume a single photon detector type for the whole ERICH (DRICH),
      // which must be a reasonable assumption; eventually may want to generalize;
      const auto pd = m_IrtDet->m_PhotonDetectors[0];
      photon->SetPhotonDetector(pd);
      {
	uint64_t vcopy = hit.cellID() & m_ReadoutCellMask;

	photon->SetVolumeCopy(vcopy);

	// Start vertex and momentum; 
	photon->SetVertexPosition(vtx);
	{
	  const auto &p = phtrack.ps();
	  photon->SetVertexMomentum(TVector3(p.x, p.y, p.z));
	}
	
	{
	  particle->FindRadiatorHistory(radiator)->AddOpticalPhoton(photon);
	  
	  // Retrieve a refractive index estimate; it is not exactly the one, which 
	  // was used in GEANT, but should be very close;
	  double ri;
	  GetFinelyBinnedTableEntry(radiator->m_ri_lookup_table, eVenergy, &ri);
	  photon->SetVertexRefractiveIndex(ri);
	}
      }

      // Hit location;
      {
	// FIXME: take sensor plane orientation into account, and also use floor()
	// rather than gaussian smearing;
	double sigma = 3.4/sqrt(12.0);
	const auto &x = hit.position();
	photon->SetDetectionPosition(TVector3(x.x + sigma*m_rngUni(), x.y + sigma*m_rngUni(), x.z));
      }     
      
      // At this point all of the CherenkovPhoton class internal variables are actually 
      // set the same way as in a standalone G4 code, except for the parent momentum at vertex;
      photon->SetDetected(true);
    } //for hit

    // FIXME: figure out where from do we have muons in the .hepmc files;
    if (!useful_photons_found) continue;

    // Loop through the radiators one by one, using the same set of photons;
    for(auto radiator: m_SelectedRadiators) {
      auto history = particle->FindRadiatorHistory(radiator);

      history->CalculateSteps();

      unsigned stCount = history->StepCount();
      if (!stCount) continue;
      auto fstep = history->GetStep(0), lstep = history->GetStep(stCount-1);

      // Let's say these are the best guesses;
      unsigned ifsec = m_IrtDet->GetSector(fstep->GetPosition());
      unsigned ilsec = m_IrtDet->GetSector(lstep->GetPosition());

      radiator->ResetLocations();

      // Give the algorithm radiator surface boundaries as encoded in D(E)Rich_geo.cpp;  
      auto s1 = radiator->GetFrontSide(ifsec);
      auto s2 = radiator->GetRearSide (ilsec);

      TVector3 from, to, p0;
      
      // FIXME: need it not at vertex, but in the radiator; as coded here, this can 
      // hardly work once the magnetic field is turned on; but this is a best guess 
      // if nothing else is available;
      const auto &p = mctrack.ps();
      
      p0 = TVector3(p.x, p.y, p.z);
      
      auto x0 = fstep->GetPosition();
      auto n0 = fstep->GetDirection(); 
      // Go backwards and ignore surface orientation mismatch;
      bool b1 = s1->GetCrossing(x0, -1*n0, &from, false);

      auto x1 = lstep->GetPosition();
      auto n1 = lstep->GetDirection();
      bool b2 = s2->GetCrossing(x1,    n1, &to);

      if (b1 && b2) {
	unsigned zbins = radiator->GetTrajectoryBinCount();
	TVector3 nn = (to - from).Unit();
	from += (0.010)*nn;
	to   -= (0.010)*nn;
	
	// FIXME: assume a straight line to the moment;
	auto span = to - from;
	double tlen = span.Mag();
	double step = tlen / zbins;
	for(unsigned iq=0; iq<zbins+1; iq++) 
	  radiator->AddLocation(from + iq*step*nn, p0);
      } //if
    } //for radiator

      // Now that all internal mctrack-level structures are populated, run IRT code;
    {
      CherenkovPID pid;
      // Should suffice for now: e+/pi+/K+/p?; FIXME: hardcoded;
      int pdg_table[] = {-11, 211, 321, 2212};
      for(unsigned ip=0; ip<sizeof(pdg_table)/sizeof(pdg_table[0]); ip++) {
	const auto &service = m_pidSvc->particle(pdg_table[ip]);
	
	pid.AddMassHypothesis(service.mass);
      } //for ip
      
	// Eventually call the IRT code;
      particle->PIDReconstruction(pid);
      
      // Populate the output collection;
      auto cbuffer = cpid.create();

      for(unsigned ir=0; ir< m_SelectedRadiators.size(); ir++) {
	const auto radiator = m_SelectedRadiators[ir];
	
	for(unsigned ip=0; ip<sizeof(pdg_table)/sizeof(pdg_table[0]); ip++) {
	  auto hypo = pid.GetHypothesis(ip);
	  eic::CherenkovPdgHypothesis hypothesis;
	  
	  // Flip sign if needed; FIXME: use reconstructed one?;
	  if (m_pidSvc->particle(mctrack.pdgID()).charge < 0.0) pdg_table[ip] *= -1;
	  
	  // Storage model is not exactly efficient, but it is simple for users;
	  hypothesis.radiator = ir;
	  hypothesis.pdg      = pdg_table[ip];
	  hypothesis.npe      = hypo->GetNpe   (radiator);
	  hypothesis.weight   = hypo->GetWeight(radiator);
	  
	  cbuffer.addoptions(hypothesis);
	} //for ip

	// Theta angle estimates, per radiator;
	{
	  unsigned npe = 0;
	  double theta = 0.0;
	  double ri = 0.0;
	  eic::CherenkovThetaAngleMeasurement rdata;

	  rdata.radiator = ir;
	  
	  auto history = particle->FindRadiatorHistory(radiator);
	  // This loop goes only over the photons which belong to a given radiator;
	  for(auto photon: history->Photons()) {
	    bool selected = false;
	    // Check whether this photon was selected by at least one mass hypothesis;
	    for(auto entry: photon->_m_Selected)
	      if (entry.second == radiator) {
		selected = true;
		break;
	      } //if

	    if (selected) {
	      npe++;
	      theta += photon->_m_PDF[radiator].GetAverage();
	      ri    += photon->GetVertexRefractiveIndex();
	    } //if
	  } //for photon

	  rdata.npe    = npe;
	  rdata.theta  = npe ? theta / npe : 0.0;
	  rdata.rindex = npe ? ri    / npe : 0.0;
	  //printf("@@@ %7.2f\n", 1000*rdata.theta);
	  cbuffer.addangles(rdata);
	}
      } //for ir

      // FIXME: well, and what does go here instead of 0?;
      cbuffer.ID({0, algorithmID()});
	
      // Reference to either MC track or a reconstructed track;
#ifdef _USE_RECONSTRUCTED_TRACKS_
      cbuffer.recID(rctrack.ID());
#else
      cbuffer.recID(mctrack.ID());
#endif
    }  
  } //for rctrack (mctrack)

  event->Reset();
  delete event;

  return StatusCode::SUCCESS;
} // Jug::PID::IRTAlgorithm::execute()

// -------------------------------------------------------------------------------------

StatusCode Jug::PID::IRTAlgorithm::finalize( void ) 
{
  info() << "IRTAlgorithm: Finalizing..." << endmsg;

  return Algorithm::finalize(); // must be executed last
} // Jug::PID::IRTAlgorithm::finalize()

// -------------------------------------------------------------------------------------


