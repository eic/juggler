
#include <map>
#include <fmt/core.h>

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
  : GaudiAlgorithm(name, svcLoc), m_IrtGeo(0), m_IrtDet(0), 
    m_ReadoutCellMask(0x0)
{
  declareProperty("inputMCParticles",                 m_inputMCParticles,              "MCParticles");
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
      std::make_unique<DataHandle<edm4hep::SimTrackerHitCollection>>((dname + "Hits").c_str(), 
									 Gaudi::DataHandle::Reader, this);
    // Output PID info collection;
    m_outputCherenkovPID =
      std::make_unique<DataHandle<edm4eic::CherenkovParticleIDCollection>>((dname + "PID").c_str(), 
								       Gaudi::DataHandle::Writer, this);
  }

    // FIXME: how about Gaudi::Property<std::vector<std::tuple>>?;
  {
    const auto &radiators = m_RadiatorConfig.value();
    
    for(auto rptr: radiators) {
      // parse radiator specification string; FIXME: expected format could be more consistent, e.g. each token is `this=that`
      //                                      FIXME: check for incorrect formatting
      /* expected tokens (delimited by space, order does not matter):
       * - "%s"          parses to  `name`
       * - "zbins=%u"    parses to  `zbins`
       * - "smearing=%s" parses to  `mode`
       * - "%fmrad"      parses to  `mrad`
       */

      // radiator specifications, to be parsed from options file
      std::string name = "";
      std::string mode = "";
      unsigned zbins = 0;
      float qmrad = -1.;
      float ri = -1.;
      float attenuation = -1.;

      // parse
      std::stringstream rptrStream(rptr);
      std::string rptrTok;
      std::regex rmEq("^.*=");
      while(rptrStream >> rptrTok) {
	if(rptrTok.find("zbins")!=std::string::npos) 
	  zbins = std::stoul(std::regex_replace(rptrTok,rmEq,""));
	else if(rptrTok.find("rindex")!=std::string::npos) 
	  ri = std::stof(std::regex_replace(rptrTok,rmEq,""));
	else if(rptrTok.find("attenuation[mm]")!=std::string::npos) 
	  attenuation = std::stof(std::regex_replace(rptrTok,rmEq,""));
	else if(rptrTok.find("smearing")!=std::string::npos) 
	  mode = std::regex_replace(rptrTok,rmEq,"");
	else if(rptrTok.find("mrad")!=std::string::npos) 
	  qmrad = std::stof(std::regex_replace(rptrTok,std::regex("mrad$"),""));
	else name = rptrTok; // otherwise assume it's the name
      }
      //printf("@@@ %s %s %d %7.1f\n", name.c_str(), mode.c_str(), zbins, qmrad);
      
      if (!zbins || qmrad<0.0 || ri<0.0 || /*attenuation<0.0 ||*/ name.compare("")==0 || mode.compare("")==0) {
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
	if (attenuation > 0.0)
	  radiator->SetReferenceAttenuationLength(attenuation);

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
  std::map<eic::Index, const edm4eic::ConstTrajectory*> rc2trajectory;
  for(const auto &trajectory: trajectories)
    rc2trajectory[trajectory.trackID()] = &trajectory;
#endif
#ifdef _USE_ON_THE_FLY_TRAJECTORIES_
  const TrajectoriesContainer* trajectories = m_inputTrajectories.get();
#endif
  // Output collection(s);
  auto &cpid               = *m_outputCherenkovPID->createAndPut();
    
  // An interface variable; FIXME: check memory cleanup;
  auto event = new CherenkovEvent(); // TODO

  //printf("%3ld track(s) and %4ld hit(s)\n", mctracks.size(), hits.size());

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
    //printf("AT least MC Loop!\n");
    // At this point have a reference to a 'mctrack', either this or that way;
    // Now just follow the logic of TrackParamTruthInit.cpp; 

    // generatorStatus = 1 means thrown G4Primary, but dd4gun uses generatorStatus == 0
    //if(mctrack.getPDG()>0) printf("@@@@ PDG: %d MomZ : %7.2f; GenStat %d, size: %d\n", mctrack.getPDG(), (mctrack.getMomentum()).z ,mctrack.getGeneratorStatus(),mctrack.parents_size());
    if (mctrack.getGeneratorStatus() != 1 ) {
      //printf("@@@@ PDG: %d MomZ : %7.2f; GenStat %d\n", mctrack.getPDG(), (mctrack.getMomentum()).z ,mctrack.getGeneratorStatus());
      if (msgLevel(MSG::DEBUG)) {
	debug() << "ignoring particle with generatorStatus = " << mctrack.getGeneratorStatus() << endmsg;
      }
      continue;
    } //if

    const double charge = m_pidSvc->particle(mctrack.getPDG()).charge;
    if (abs(charge) < std::numeric_limits<double>::epsilon()) {
      if (msgLevel(MSG::DEBUG)) {
	debug() << "ignoring neutral particle" << endmsg;
      }
      continue;
    } //if

      // FIXME: consider only primaries for the time being;
    //if (mctrack.parents_size()!=2) continue;
    //printf("@@@  %3d %d\n", mctrack.getPDG(), mctrack.parents_size());
    // FIXME: a hack; remove muons;
    if (mctrack.getPDG() == 13) continue;
    
    //printf("CAME HERE  MC Loop!\n");
    auto particle = new ChargedParticle(mctrack.getPDG()); // TODO
    //printf("HepMC Particle PDG %d \n",mctrack.getPDG());
    event->AddChargedParticle(particle);  // TODO

    // Start radiator history for all known radiators;
    for(auto rptr: m_IrtDet->Radiators())
      particle->StartRadiatorHistory(std::make_pair(rptr.second, new RadiatorHistory())); // TODO
      
    // Distribute photons over radiators where they were presumably produced, to 
    // create a structure which would mimic a standalone G4 stepper behavior;
    bool useful_photons_found = false;
    //printf("%3d\n", mctracks.size());
    for(const auto &hit: hits) {
      // FIXME: range checks;
      const auto &phtrack = hit.getMCParticle();
      if(phtrack.getPDG()!=-22) {
        fprintf(stderr,"ERROR: this hit was not from an optical photon: PDG = %d\n",phtrack.getPDG());
        continue;
      }

      // FIXME: yes, use MC truth here; not really needed I guess;
      // FIXME: range checks;
      //printf("Here: %d %d %3d!\n", phtrack.parents()[0], mctrack.ID(), hit.truth().trackID);
      //if (phtrack.parents()[0] != mctrack.ID()) continue;
      
      // Vertex where photon was created;
      const auto &vs = phtrack.getVertex();
      TVector3 vtx(vs.x, vs.y, vs.z);
      // printf("@@@ photon vtx: %f %f %f\n", vs.x, vs.y, vs.z);
      {
        // const auto &ps = phtrack.getMomentum();
        // double phi = TVector3(ps.x, ps.y, ps.z).Phi();
        // double theta = TVector3(ps.x, ps.y, ps.z).Theta();
        // printf("@@@ theta, phi = %7.2f %7.2f\n",theta,phi);
        //if (fabs(phi - M_PI/2) > M_PI/4 && fabs(phi + M_PI/2) > M_PI/4) continue;
      }
      //if (vs.z < 2500 || vs.z > 2700) continue;
      
      // FIXME: this is in general correct, but needs refinement;
      TVector3 ip(0,0,0);
      // FIXME: unify with the code below;
      auto radiator = m_IrtDet->GuessRadiator(vtx, (vtx - ip).Unit());
      // FIXME: do it better later;
      if (!radiator) continue;
      //printf("Here-1\n");

      // Simulate QE & geometric sensor efficiency; FIXME: hit.energy() is numerically 
      // in GeV units, but Gaudi::Units::GeV = 1000; prefer to convert photon energies 
      // to [eV] in all places by hand;
      double eVenergy = 1E9*hit.getEDep();
      //printf("%f\n", eVenergy);
      if (!QE_pass(eVenergy, m_rngUni()) || 
	  m_rngUni() > /*m_GeometricEfficiency.value()**/m_SafetyFactor.value()) {

	// Prefer to create the stepping history in a clean way, namely use only those 
	// photons which were rejected; FIXME: optionally include all; anyway, now 
	// figure out (1) to which detector sector does this 3D point belongs, (2) 
	// which radiator this could have been;

	// Add the point to the radiator history buffer; FIXME: all this is not 
	// needed if ACTS trajectory parameterization is used;
	particle->FindRadiatorHistory(radiator)->AddStepBufferPoint(phtrack.getTime(), vtx); // TODO? maybe not this way, this is photon pinning
                                                                                             // TODO use TrackSegment's TrackPoints instead

	//printf("Here-2A\n");
	continue;
      } //if
      //printf("Here-2B\n");

      //if (vs.z < 2500 || vs.z > 2700) continue;
      //if (vs.z < 2200 || vs.z > 2400) continue;

      // Hit location after pixelization;
      // FIXME: use cellID decoder, but we need gaps... see https://github.com/AIDASoft/DD4hep/issues/958
      TVector3 lxyPixel;
      // FIXME: '0' - for now assume a single photon detector type for the whole ERICH (DRICH),
      // which must be a reasonable assumption; eventually may want to generalize;
      const auto pd = m_IrtDet->m_PhotonDetectors[0]; // TODO
      uint64_t vcopy = hit.getCellID() & m_ReadoutCellMask; // TODO
      {
	const auto irt = pd->GetIRT(vcopy); // TODO
	auto sensor = dynamic_cast<const FlatSurface*>(irt->tail()->GetSurface()); // TODO

	double pitch = m_SensorPixelPitch.value(), sens = m_SensorPixelSize.value();//pitch - gap;
	const auto &x = hit.getPosition();
	TVector3 lpt(x.x, x.y, x.z);
	double lxy[2] = {sensor->GetLocalX(lpt), sensor->GetLocalY(lpt)}, lxyPixels[2] = {0.0};
	//printf("%7.2f %7.2f\n", lxy[0], lxy[1]);
	int ixy[2];
	for(unsigned iq=0; iq<2; iq++) {
	  ixy[iq] = (int)floor(lxy[iq]/pitch);//m_SensorPixelPitch.value());
	  
	  // "+0.5" assumes even number of pixels in each projection;
	  lxyPixels[iq] = pitch*(ixy[iq] + 0.5);
	} //for iq

	// Well, just ignore this small fraction of unlucky photons (do not even use them for track pinning);
	if (fabs(lxy[0] - lxyPixels[0]) > sens/2 || fabs(lxy[1] - lxyPixels[1]) > sens/2) continue;

	lxyPixel = sensor->GetSpacePoint(lxyPixels[0], lxyPixels[1]); // TODO: what is this number? do we need PhotoMultiplierHitDigi's pixel hit position?
      }     

      useful_photons_found = true;

      // Photon accepted; add it to the internal event structure; // TODO
      auto photon = new OpticalPhoton();
      photon->SetVolumeCopy(vcopy);
      photon->SetDetectionPosition(lxyPixel);
      photon->SetPhotonDetector(pd);

      // EICrecon cross check:
      fmt::print("[cross-check] cell_id={:#X}  copy={}  pixel_pos=( {:>10.2f} {:>10.2f} {:>10.2f} )\n",
          hit.getCellID(), vcopy, lxyPixel.x(), lxyPixel.y(), lxyPixel.z());

      {
	// Start vertex and momentum; 
	photon->SetVertexPosition(vtx);
	{
	  const auto &p = phtrack.getMomentum();
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
      
      // At this point all of the CherenkovPhoton class internal variables are actually 
      // set the same way as in a standalone G4 code, except for the parent momentum at vertex;
      photon->SetDetected(true);
    } //for hit

    // FIXME: figure out where from do we have muons in the .hepmc files;
    if (!useful_photons_found) continue;

    // Loop through the radiators one by one;
    for(auto radiator: m_SelectedRadiators) {
      auto history = particle->FindRadiatorHistory(radiator);

      history->CalculateSteps(); // TODO: NOT NECESSARY ALREADY DONE

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
      const auto &p = mctrack.getMomentum();
      
      p0 = TVector3(p.x, p.y, p.z);
      
      auto x0 = fstep->GetPosition();
      auto n0 = fstep->GetDirection(); 
      // Go backwards and ignore surface orientation mismatch;
      bool b1 = s1->GetCrossing(x0, -1*n0, &from, false);

      auto x1 = lstep->GetPosition();
      auto n1 = lstep->GetDirection();
      bool b2 = s2->GetCrossing(x1,    n1, &to);

      if (b1 && b2) {
	const unsigned zbins = radiator->GetTrajectoryBinCount();
	TVector3 nn = (to - from).Unit();
	// FIXME: hardcoded;
	from += (0.010)*nn;
	to   -= (0.010)*nn;
	
#if 1
	// No assumptions about straight or curved trajectory;
	{
	  double zfrom = from.z(), zto = to.z();
	  //double step = fabs(zto - zfrom) / zbins;
	  double step = (zto - zfrom) / zbins;

	  for(unsigned iq=0; iq<zbins+1; iq++) {
	    double z = zfrom + iq*step;

	    //printf("%d %f\n", iq, z);
	    // Find two neighboring points in the history array, which would 
	    // allow one to compute the coordinate at this z in a best way;
	    // FIXME: this loop is not dramatically efficient;
	    for(unsigned is=1; is<stCount; is++) {
	      auto st1 = history->GetStep(is-1), st2 = history->GetStep(is);
	      const auto &pt1 = st1->GetPosition(), &pt2 = st2->GetPosition();
	      double z1 = pt1.z(), z2 = pt2.z();

	      //printf("%3d -> %f %f\n", is, z1, z2);
	      if (fabs(z) < fabs(z1) || (fabs(z) >= fabs(z1) && fabs(z) < fabs(z2)) || is == stCount-1) {
		auto p12 = pt2 - pt1, n12 = p12.Unit();
		//double len = (pt2 - pt1).Mag();
		//double zlen = fabs((pt2 - pt1).z());//Mag();
		double zlen = (pt2 - pt1).z();//Mag();

		//radiator->AddLocation(pt1 + (fabs(z - z1)/len)*nn, p0);//.Mag()*n12);
		radiator->AddLocation(pt1 + ((z - z1)/zlen)*p12, p0.Mag()*n12); // TODO: can we just use the TrackPoint position?
		//printf("(1) %f %f %f\n", 
		//     (pt1 + ((z - z1)/zlen)*p12).x(), 
		//     (pt1 + ((z - z1)/zlen)*p12).y(), 
		//     (pt1 + ((z - z1)/zlen)*p12).z());
		break;
	      } //if
	    } //for is
	  } //for iq
	}
#else
	{
	  // FIXME: assume a straight line to the moment;
	  auto span = to - from;
	  double tlen = span.Mag();
	  double step = tlen / zbins;
	  for(unsigned iq=0; iq<zbins+1; iq++) {
	    radiator->AddLocation(from + iq*step*nn, p0);
	    printf("(2) %f %f %f\n", 
		   (from + iq*step*nn).x(),
		   (from + iq*step*nn).y(),
		   (from + iq*step*nn).z());
	  } //for iq
	}
#endif
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
	  edm4eic::CherenkovPdgHypothesis hypothesis;
	  //printf("HEPMC Entering ip loop \n"); 
	  // Flip sign if needed; FIXME: use reconstructed one?;
	  if (m_pidSvc->particle(mctrack.getPDG()).charge < 0.0) pdg_table[ip] *= -1;
	  
	  // Storage model is not exactly efficient, but it is simple for users;
	  hypothesis.radiator = ir;
	  hypothesis.pdg      = pdg_table[ip];
	  hypothesis.npe      = hypo->GetNpe   (radiator);
	  hypothesis.weight   = hypo->GetWeight(radiator);
	  
	  cbuffer.addToOptions(hypothesis);

	} //for ip

	// Theta angle estimates, per radiator;
	{
	  unsigned npe = 0;
	  double theta = 0.0;
	  double ri = 0.0;
	  double wl = 0.0;
	  double wtsum = 0.0;
          std::vector<double> thphot;
          std::vector<double> phiphot;
	  edm4eic::CherenkovThetaAngleMeasurement rdata;

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

	     phiphot.push_back(photon->m_Phi[radiator]);
	    //printf("PHI ANGLE: ###  %f\n", phi);
	    if (selected) {
	      //if (selected &&
	      //(fabs(phi - M_PI/2) < M_PI/4 || fabs(phi + M_PI/2) < M_PI/4)) {
	      npe++;
	      double wt = 1.0;//fabs(sin(photon->m_Phi[radiator]));
	      wtsum += wt;
	      theta += wt*photon->_m_PDF[radiator].GetAverage();
              thphot.push_back(photon->_m_PDF[radiator].GetAverage());

              cbuffer.addToThetaphoton(photon->_m_PDF[radiator].GetAverage());
              cbuffer.addToPhiphoton(photon->m_Phi[radiator]);

	      ri    +=    photon->GetVertexRefractiveIndex();
	      // FIXME: hardcoded;
	      wl    +=    1239.8/(1E9*photon->GetVertexMomentum().Mag());
	    } //if
	  } //for photon
          /*double phiph=0.;
          if(phiphot.size()!= thphot.size()) printf("Phot theta phi size mismatch! \n");
          for(unsigned int i=0; i<thphot.size(); i++){
            if(57.29*phiphot[i]<0)  phiph = 360+ 57.29*phiphot[i];
             else phiph = 57.29*phiphot[i];
            if(wtsum>0) printf("@@ Npe: %d Ring Theta: %7.2f %d th Photon theta %7.2f phi: %7.2f \n",npe, 1000*theta / wtsum,i,1000*thphot[i], phiph  );

          }
          */
	  rdata.npe        = npe;
	  //rdata.theta  = npe ? theta / npe : 0.0;
	  rdata.theta      = wtsum ? theta / wtsum : 0.0;
	  rdata.rindex     = npe   ? ri    / npe   : 0.0;
	  rdata.wavelength = npe   ? wl    / npe   : 0.0;
          //rdata.thetaphoton = thphot;
          //rdata.phiphoton = phiphot;
	  //printf("PHOTON ANGLE : @@@ %7.2f\n", 1000*rdata.theta);
	  cbuffer.addToAngles(rdata);
	}
      } //for ir

      // FIXME: well, and what does go here instead of 0?;
      // cbuffer.ID({0, algorithmID()});
	
      // Reference to either MC track or a reconstructed track;
#ifdef _USE_RECONSTRUCTED_TRACKS_
      cbuffer.recID(rctrack.ID()); // do not use, out of date
#else
      cbuffer.setAssociatedParticle(mctrack);
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


