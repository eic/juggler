
#include <map>

#include "TFile.h"

#include "DDRec/CellIDPositionConverter.h"

#include "IRTAlgorithm.h"

#include "IRT/CherenkovEvent.h"
#include "IRT/CherenkovDetectorCollection.h"

// -------------------------------------------------------------------------------------

Jug::Reco::IRTAlgorithm::IRTAlgorithm(const std::string& name, ISvcLocator* svcLoc) 
  : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()), m_IrtGeo(0), m_IrtDet(0) 
{
  declareProperty("inputMCParticles",                 m_inputMCParticles,              "");
  declareProperty("inputHitCollection",               m_inputHitCollection,            "");
#ifdef _USE_RECONSTRUCTED_TRACKS_
  declareProperty("inputRecoParticles",               m_inputRecoParticles,            "");
#endif
  declareProperty("inputTrajectories",                m_inputTrajectories,             "");

  declareProperty("CherenkovPID",                     m_outputCherenkovPID,            "");
} // Jug::Reco::IRTAlgorithm::IRTAlgorithm()

// -------------------------------------------------------------------------------------

StatusCode Jug::Reco::IRTAlgorithm::initialize( void ) 
{
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
    auto const &config = m_ConfigFile.value();

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

      m_IrtDet = m_IrtGeo->GetDetector(0);
    } else {
      // Otherwise create IRT detector collection geometry;
      m_IrtGeo = new CherenkovDetectorCollection();
      m_IrtGeo->AddNewDetector();

      m_IrtDet = m_IrtGeo->GetDetector(0);
      // TODO: copy over Alexander's FIXME comments
      
      // get detector elements
      auto dd4Det = m_geoSvc->detector();
      auto erichDet = dd4Det->detector("ERICH");
      //auto drichDet = dd4Det->detector("DRICH");
      
      // get detector positions
      auto erichPos = erichDet.placement().position();
      //auto drichPos = drichDet.placement().position();
      
      // set IRT container volume
      auto boundary = new FlatSurface(TVector3(0,0,0), TVector3(1,0,0), TVector3(0,1,0));
      m_IrtGeo->SetContainerVolume(m_IrtDet, (G4LogicalVolume*)(0x0), 0, boundary);
      
      // eRICH loop :::::::::::::::::::::::
      int sector;
      for(auto const& [erichDEname, erichDE] : erichDet.children()) {
	//info() << "FOUND ERICH DE: " << erichDEname << endmsg;
	auto pos = erichPos + erichDE.placement().position();
	
	// aerogel and filter
	if(erichDEname.find("aerogel")!=std::string::npos || erichDEname.find("filter")!=std::string::npos) {
	  double thickness = 2 * erichDE.volume().boundingBox().dimensions()[2];
	  sector = erichDE.id();
	  info() << "ERICH RADIATOR: " << erichDEname
		 << "\n\t(x,y,z)-position = " << pos.x() << ", " << pos.y() << ", " << pos.z()
		 << "\n\tsector = " << sector
		 << "\n\tthickness = " << thickness
		 << endmsg;
	  if(sector==0) {
	    if(erichDEname.find("aerogel")!=std::string::npos) {
	      auto aerogelSurf = new FlatSurface(
						 (1/mm)*TVector3(0,0,pos.z()), // TODO: correct position?
						 TVector3(1,0,0),
						 TVector3(0,1,0)
						 );
	      m_IrtGeo->AddFlatRadiator(m_IrtDet, (G4LogicalVolume*)(0x1), 0, aerogelSurf, thickness/mm); // TODO: correct units?
	    } else {
	      auto filterSurf = new FlatSurface(
						(1/mm)*TVector3(0,0,pos.z()-0.01), // TODO: correct position?
						TVector3(1,0,0),
						TVector3(0,1,0)
						);
	      m_IrtGeo->AddFlatRadiator(m_IrtDet, (G4LogicalVolume*)(0x2), 0, filterSurf, thickness/mm);
	    }
	  }
	} // end if aerogel
	
	
	// sensors
	if(erichDEname.find("sensor")!=std::string::npos) {
	  info() << "ERICH SENSOR: " << erichDEname
		 << "\n\t(x,y,z)-position = " << pos.x() << ", " << pos.y() << ", " << pos.z()
		 << "\n\tid = " << erichDE.id()
		 << endmsg;
	  auto sensorSurf = new FlatSurface(
					    (1/mm)*TVector3(0.0, 0.0, pos.z()), // TODO: why not `pos.x(), pos.y(), pos.z()`?
					    TVector3(1,0,0),
					    TVector3(0,1,0)
					    );
	  m_IrtDet->AddPhotonDetector(new CherenkovPhotonDetector(0,0,sensorSurf));
	}
	
      } // end loop over eRICH detector elements
    } //if
  }

#if _TODAY_  
  // set radiator refractive indices
  info() << "C4F10 RINDEX: " << endmsg;
  auto gasRmatrix = dd4Det->material("C4F10_ERICH").property("RINDEX"); // TODO: get material from gas volume instead
  for(size_t r=0; r<gasRmatrix->GetRows(); r++) {
    info() << "   " << gasRmatrix->Get(r,0) << "   " << gasRmatrix->Get(r,1) << endmsg;
  };
#endif
  // TODO: irtDet->Radiators()[0]->SetReferenceRefractiveIndex( N ) // aerogel
  // TODO: irtDet->Radiators()[1]->SetReferenceRefractiveIndex( N ) // filter
  // TODO: irtDet->Radiators()[2]->SetReferenceRefractiveIndex( N ) // gas

  // Need a random number generator for the QE stuff;
  {
    auto randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
    auto sc = m_rngUni.initialize(randSvc, Rndm::Flat(0., 1.));
    if (!sc.isSuccess()) {
      error() << "Cannot initialize random generator!" << endmsg;
      return StatusCode::FAILURE;
    } //if
  }

  // Initialize the QE lookup table;
  configure_QE_lookup_table(m_QE_input_data.value(), m_QEbins.value());

  return StatusCode::SUCCESS;
} // Jug::Reco::IRTAlgorithm::initialize()

// -------------------------------------------------------------------------------------

StatusCode Jug::Reco::IRTAlgorithm::execute( void )
{
  // Input collection(s);
  const auto &hits         = *m_inputHitCollection.get();
  const auto &mctracks     = *m_inputMCParticles.get();
#ifdef _USE_RECONSTRUCTED_TRACKS_
  const auto &rctracks     = *m_inputRecoParticles.get();
#endif
  const auto &trajectories = *m_inputTrajectories.get();

  // Output collection(s);
  auto &cpid          = *m_outputCherenkovPID.createAndPut();

  // First populate the trajectory-to-reconstructed (or -to-simulated) mapping table;
  std::map<eic::Index, const eic::ConstTrajectory*> rc2trajectory;
  for(const auto &trajectory: trajectories)
    rc2trajectory[trajectory.trackID()] = &trajectory;
    
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
    // Follow the logic of TrackParamTruthInit.cpp; 

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

    auto particle = new ChargedParticle(mctrack.pdgID());
    event->AddChargedParticle(particle);

    // FIXME: hardcoded;
    auto aerogel = m_IrtDet->Radiators()[1];
    aerogel->m_AverageRefractiveIndex = aerogel->n();

    {
      auto history = new RadiatorHistory();
      // FIXME: yes, for now assume all the photons were produced in aerogel; 
      particle->StartRadiatorHistory(std::make_pair(aerogel, history));

      {
	// FIXME: hardcoded; give the algorithm aerogel surface boundaries as encoded in ERich_geo.cpp;  
	auto s1 = m_IrtDet->m_OpticalBoundaries[1], s2 = m_IrtDet->m_OpticalBoundaries[2];
	//printf("%ld\n", m_IrtDet->m_OpticalBoundaries.size());

	// If trajectory parameterizations are actually available, use them;
#ifdef _USE_RECONSTRUCTED_TRACKS_
	auto index = rctrack.ID();
#else
	auto index = mctrack.ID();
#endif
	auto trajectory = rc2trajectory.find(index) == rc2trajectory.end() ? 0 : 
	  rc2trajectory[index];

	// FIXME: merge this if-elseif somehow; 
	if (trajectory && trajectory->points().size()) {
	  // Assume 's1' (front surface) is good enough;
	  auto surface = s1->GetSurface();
	  const eic::TrajectoryPoint *best = 0;

	  // In case of aerogel, just find the closest point; gas case will be more complicated,
	  // but in general pretty much straightforward too;
	  for(const auto &point: trajectory->points()) 
	    // FIXME: optimize for efficiency later; now calculate every time new;
	    if (!best || IRTAlgorithmServices::GetDistance(surface, &point) < 
		IRTAlgorithmServices::GetDistance(surface, best))
	      best = &point;

	  // FIXME: well, may want to check how far away the best point actually is;
	  {
	    TVector3 from, to;
	    IRTAlgorithmServices::GetCrossing(s1->GetSurface(), best, &from);
	    IRTAlgorithmServices::GetCrossing(s2->GetSurface(), best, &to);
	    
	    auto p0 = IRTAlgorithmServices::GetMomentum(best);
	    history->AddStep(new ChargedParticleStep(from, p0));
	    history->AddStep(new ChargedParticleStep(  to, p0));
	  }
	} else {
	  // FIXME: need it not at vertex, but in the radiator; as coded here, this can 
	  // hardly work once the magnetic field is turned on; but this is a best guess 
	  // if nothing else is available;
	  const auto &vtx = mctrack.vs(), &p = mctrack.ps();
  
	  auto x0 = TVector3(vtx.x, vtx.y, vtx.z), p0 = TVector3(p.x, p.y, p.z), n0 = p0.Unit();
		
	  TVector3 from, to;
	  s1->GetSurface()->GetCrossing(x0, n0, &from);
	  s2->GetSurface()->GetCrossing(x0, n0, &to);
	  
	  history->AddStep(new ChargedParticleStep(from, p0));
	  history->AddStep(new ChargedParticleStep(  to, p0));
	} //if
      }

      for(const auto &hit: hits) {
	// FIXME: yes, use MC truth here; not really needed I guess; 
	if (hit.g4ID() != mctrack.ID()) continue;

	// Simulate QE; FIXME: random seed, reproducibility, etc;
	// FIXME: check units;
	if (!QE_pass(hit.energy(), m_rngUni()*m_GeometricEfficiency.value())) continue;
	
	// Photon accepted; add it to the internal event structure;
	auto photon = new OpticalPhoton();
	  
	{
	  const auto &x = hit.position();
	  photon->SetDetectionPosition(TVector3(x.x, x.y, x.z));
	}

	//
	// FIXME: detector ID check needed?;
	// 

	//<id>system:8,module:12,x:20:16,y:16</id>
	//printf("%4ld %4ld %4ld %4ld\n", 
	//     ( hit.cellID        & 0x00FF), 
	//     ((hit.cellID >>  8) & 0x0FFF),
	//     ((hit.cellID >> 12) & 0xFFFF),
	//     ((hit.cellID >> 28) & 0xFFFF));
	unsigned module = (hit.cellID() >>  8) & 0x0FFF;
	if (module < m_IrtDet->m_PhotonDetectors.size()) {
	  photon->SetPhotonDetector(m_IrtDet->m_PhotonDetectors[module]);
	  photon->SetDetected(true);
	  
	  history->AddOpticalPhoton(photon);
	} //if
      } //for hit

	//printf("%d %d %d\n", mctrack.ID, mctrack.pdgID, mctrack.g4Parent);

      // Now that all internal mctrack-level structures are populated, run IRT code;
      {
	// IRT side of the story;
	CherenkovPID pid;
	// Should suffice for now: e+/pi+/K+/p?;
	int pdg_table[] = {-11, 211, 321, 2212};
	// Flip sign if needed; FIXME: use reconstructed one?;
	for(unsigned ip=0; ip<sizeof(pdg_table)/sizeof(pdg_table[0]); ip++) {
	  const auto &service = m_pidSvc->particle(mctrack.pdgID());
	  if (service.charge < 0.0) pdg_table[ip] *= -1;

	  pid.AddMassHypothesis(service.mass);
	} //for ip

	// Eventually call the IRT code;
	particle->PIDReconstruction(pid);

	// Now the interesting part comes: how to populate the output collection;
	{
	  auto cbuffer = cpid.create();

	  for(unsigned ip=0; ip<sizeof(pdg_table)/sizeof(pdg_table[0]); ip++) {
	    auto hypo = pid.GetHypothesis(ip);
	    eic::CherenkovPdgHypothesis hypothesis;

	    hypothesis.pdg    = pdg_table[ip];
	    hypothesis.npe    = hypo->GetNph();
	    hypothesis.weight = hypo->GetWeight();

	    cbuffer.addoptions(hypothesis);
	  } //for ip

	  // FIXME: well, and what does go here instead of 0?;
	  cbuffer.ID({0, algorithmID()});

	  // Reference to either MC track or a reconstructed track;
#ifdef _USE_RECONSTRUCTED_TRACKS_
	  cbuffer.recID(rctrack.ID());
#else
	  cbuffer.recID(mctrack.ID());
#endif

	  //cpid.push_back(cbuffer);
	}
      }
    }  
  } //for rctrack (mctrack)

  event->Reset();
  delete event;

  return StatusCode::SUCCESS;
} // Jug::Reco::IRTAlgorithm::execute()

// -------------------------------------------------------------------------------------

StatusCode Jug::Reco::IRTAlgorithm::finalize( void ) 
{
  info() << "IRTAlgorithm: Finalizing..." << endmsg;

  return Algorithm::finalize(); // must be executed last
} // Jug::Reco::IRTAlgorithm::finalize()

// -------------------------------------------------------------------------------------
