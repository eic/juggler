
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
#ifdef _USE_TRAJECTORIES_
  declareProperty("inputTrajectories",                m_inputTrajectories,             "");
#endif
} // Jug::PID::IRTAlgorithm::IRTAlgorithm()

// -------------------------------------------------------------------------------------

StatusCode Jug::PID::IRTAlgorithm::initialize( void ) 
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

      if (dname.empty()) {
        error() << "No RICH detector name provided." << endmsg;
        return StatusCode::FAILURE;
      } //if

    } else { // USE GEOMETRY SERVICE to help create IRT detector collection geometry

      bool debug_geosvc = 1; // print out geometry attributes, etc.

      // set IRT geometry
      m_IrtGeo = new CherenkovDetectorCollection();
      m_IrtGeo->AddNewDetector(dname.c_str());
      m_IrtDet = m_IrtGeo->GetDetector(dname.c_str());
      
      // get RICH detector element and position
      auto dd4Det = m_geoSvc->detector();
      auto richDet = dd4Det->detector(dname.c_str());
      auto richPos = richDet.placement().position();

      // RICH-specific settings
      int zDirection = (dname.compare("DRICH")==0) ? 1 : -1;
      if(debug_geosvc)
        info() << "Initialize IRT algorithm for " << (zDirection<0?"-":"+") << "z endcap RICH (" << dname << ")" << endmsg;
      int nSectors;
      TVector3 nx,ny;
      std::function<std::pair<int,int>(int)> id2secmod; // DetElement id -> pair(sector#,module#)
      switch(zDirection) {
        case -1:
          nSectors = 1;
          nx = TVector3(1,0,0);
          ny = TVector3(0,1,0);
          id2secmod = [](int id){ return std::pair<int,int>(0,id); };
          break;
        case  1:
          nSectors = 6;
          nx = TVector3(1,0,0);
          ny = TVector3(0,-1,0);
          id2secmod = [](int id){ return std::pair<int,int>(id&0x7,id>>3); }; 
                                  // FIXME: make sure this is in sync with `athena/src/*Rich_geo.cpp`
          break;
      };
      
      // set IRT container volume
      // FIXME: Z-location does not really matter here, right?; but Z-axis orientation does;
      // FIXME: have no connection to GEANT G4LogicalVolume pointers; however all is needed 
      // is to make them unique so that std::map works internally; resort to using integers, 
      // who cares; material pointer can seemingly be '0', and the effective refractive index 
      // for all radiators will be assigned at the end by hand; FIXME: should assign it on 
      // per-photon basis, at birth, like standalone GEANT code does;
      for(int isec=0; isec<nSectors; isec++) { // FIXME: do we need a sector loop? probably not...
        m_IrtGeo->SetContainerVolume(
            m_IrtDet, "GasVolume", isec,
            (G4LogicalVolume*)(0x0), 0, new FlatSurface(TVector3(0,0,0), nx, ny)
            );
      };

      // FIXME: Get access to the readout structure decoder
      // m_IrtDet->SetReadoutCellMask( ... )

      // set IRT sensors // FIXME: '0' stands for the unknown (and irrelevant) G4LogicalVolume;
      auto pd = new CherenkovPhotonDetector(0, 0);
      m_IrtGeo->AddPhotonDetector(m_IrtDet, 0, pd);
      
      // loop over RICH detector elements
      bool radClosed = false;
      for(auto const& [richDEname, richDE] : richDet.children()) {

        // get element attributes
        //if(debug_geosvc) info() << "FOUND RICH element: " << richDEname << endmsg;
        auto dePos = richPos + richDE.placement().position();
        int deID = richDE.id();

        // aerogel and filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(richDEname.find("aerogel")!=std::string::npos || richDEname.find("filter")!=std::string::npos) {
          double thickness = 2 * richDE.volume().boundingBox().dimensions()[2];
          int isec = deID;
          if(debug_geosvc) {
            info() << "RICH RADIATOR: " << richDEname
                   << "\n\t(x,y,z)-position = " << dePos.x() << ", " << dePos.y() << ", " << dePos.z()
                   << "\n\tsector = " << isec
                   << "\n\tthickness = " << thickness
                   << endmsg;
          }
          if(isec<nSectors) { // FIXME: need sector loop?
            if(richDEname.find("aerogel")!=std::string::npos) {
              auto aerogelSurf = new FlatSurface( (1/mm)*TVector3(0,0,dePos.z()), nx, ny);
              m_IrtGeo->AddFlatRadiator(m_IrtDet, "Aerogel", isec, (G4LogicalVolume*)(0x1), 0, aerogelSurf, thickness/mm);
            } else { // elif filter
              auto filterSurf = new FlatSurface( (1/mm)*TVector3(0,0,dePos.z()), nx, ny); // NOTE: there is an airgap in geometry
              m_IrtGeo->AddFlatRadiator(m_IrtDet, "Filter", isec, (G4LogicalVolume*)(0x2), 0, filterSurf, thickness/mm);
            }
          }
        }

        // mirrors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(richDEname.find("mirror")!=std::string::npos) {
          int isec = deID;

          // spherical mirror DetElement solid is a BooleanSolid; to get the sphere attributes
          // we need to access both the primitive Sphere and its relative positioning
          dd4hep::Solid sphPrim;
          auto sphPos = dePos;
          findPrimitive("TGeoSphere",richDE.solid(),sphPrim,sphPos); // get sphere primitive
          auto sph = (dd4hep::Sphere) sphPrim;

          // for some reason, the sector z-rotation is not accounted for in `findPrimitive`, so we correct for it here:
          if(richDE.placement().matrix().IsRotAboutZ()) {
            Double_t sphPosArrLocal[3], sphPosArrMaster[3];
            sphPos.GetCoordinates(sphPosArrLocal);
            richDE.placement().matrix().LocalToMaster(sphPosArrLocal,sphPosArrMaster);
            sphPos.SetCoordinates(sphPosArrMaster);
          } else error() << "richDE.placement().matrix() is not a z-rotation; cross check mirror center coords!!!" << endmsg;

          // mirror attributes
          double mirrorRadius = sph.rMin();
          if(debug_geosvc) {
            info() << "RICH MIRROR: " << richDEname
                   << "\n\t(x,y,z)-position = " << sphPos.x() << ", " << sphPos.y() << ", " << sphPos.z()
                   << "\n\tsector = " << isec
                   << "\n\tradius = " << mirrorRadius
                   << endmsg;
          }
          auto mirrorSurf = new SphericalSurface(
              (1/mm)*TVector3(sphPos.x(),sphPos.y(),sphPos.z()),
              mirrorRadius/mm
              );
          m_IrtDet->AddOpticalBoundary(isec, new OpticalBoundary(m_IrtDet->GetContainerVolume(), mirrorSurf, false));
        }

        // sensors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(richDEname.find("sensor")!=std::string::npos) {
          int ielem = deID;
          if(debug_geosvc) {
            auto secmod = id2secmod(ielem);
            int isec = secmod.first;
            int imod = secmod.second;
            info() << "RICH SENSOR: " << richDEname
                   << "\n\t(x,y,z)-position = " << dePos.x() << ", " << dePos.y() << ", " << dePos.z()
                   << "\n\tid   = " << ielem
                   << "\n\tisec = " << isec
                   << "\n\timod = " << imod
                   << endmsg;
          }

          // FIXME: why not `dePos.x(), dePos.y(), dePos.z()`? why the orientation `nx.cross(ny)`, rather than along sphere radius?
          auto sensorSurf = new FlatSurface( (1/mm)*TVector3(0.0, 0.0, dePos.z()), nx, ny);
          m_IrtDet->CreatePhotonDetectorInstance(0, pd, ielem, sensorSurf);

          // close IRT gas radiator by hand (once) // FIXME: why don't we do this after this DetElement loop?
          if(!radClosed && zDirection==-1) { // FIXME: needs to be done for dRICH too
            m_IrtDet->GetRadiator("GasVolume")->m_Borders[0].second = dynamic_cast<ParametricSurface*>(sensorSurf);
            radClosed = true;
          }
        }

      } // end loop over RICH detector elements

    } //if

      // Input hit collection;
    m_inputHitCollection =
      std::make_unique<DataHandle<dd4pod::PhotoMultiplierHitCollection>>((dname + "Hits").c_str(), 
									 Gaudi::DataHandle::Reader, this);
    // Output PID info collection;
    m_outputCherenkovPID =
      std::make_unique<DataHandle<eic::CherenkovParticleIDCollection>>((dname + "PID").c_str(), 
								       Gaudi::DataHandle::Writer, this);

    // FIXME: how about Gaudi::Property<std::vector<std::tuple>>?;
    {
      const auto &radiators = m_RadiatorConfig.value();
      
      char name[128-1], mode[128-1];
      unsigned zbins;
      float mrad;
      // FIXME: yes, a fixed format (and fixed order!) for now;
      const char *format = "%s zbins=%u smearing=%s %fmrad";
      
      for(auto rptr: radiators) {
	sscanf(rptr.c_str(), format, &name, &zbins, &mode, &mrad);
	//printf("%s %s %d %7.1f\n", name, mode, zbins, mrad);

	if (!zbins || mrad < 0.0) {
	  error() << "Provided parameters do not pass a sanity check." << endmsg;
	  return StatusCode::FAILURE;
	} //if

	// FIXME: sanity check would not hurt;
	auto radiator = m_IrtDet->GetRadiator(name);
	if (radiator) {
	  m_SelectedRadiators.push_back(radiator);
	  
	  radiator->m_AverageRefractiveIndex = radiator->n();
	  if (mrad) {
	    // Well, want them readable in the config file, but pass in [rad] to the IRT algorithm;
	    mrad /= 1000.;

	    if (strcmp(mode, "uniform") && strcmp(mode, "gaussian")) {
	      error() << "Unknown smearing mode." << endmsg;
	      return StatusCode::FAILURE;
	    } //if

	    (!strcmp(mode, "uniform") ? radiator->SetUniformSmearing(mrad) : radiator->SetGaussianSmearing(mrad));
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
  }

#if _TODAY_  
  // set radiator refractive indices
  info() << "C4F10 RINDEX: " << endmsg;
  auto gasRmatrix = dd4Det->material("C4F10_ERICH").property("RINDEX"); // FIXME: get material from gas volume instead
  for(size_t r=0; r<gasRmatrix->GetRows(); r++) {
    info() << "   " << gasRmatrix->Get(r,0) << "   " << gasRmatrix->Get(r,1) << endmsg;
  };
#endif
  // FIXME: irtDet->Radiators()[0]->SetReferenceRefractiveIndex( N ) // aerogel
  // FIXME: irtDet->Radiators()[1]->SetReferenceRefractiveIndex( N ) // filter
  // FIXME: irtDet->Radiators()[2]->SetReferenceRefractiveIndex( N ) // gas

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
#ifdef _USE_TRAJECTORIES_
  const auto &trajectories = *m_inputTrajectories.get();
#endif
  // Output collection(s);
  auto &cpid               = *m_outputCherenkovPID->createAndPut();

  // First populate the trajectory-to-reconstructed (or -to-simulated) mapping table;
#ifdef _USE_TRAJECTORIES_
  std::map<eic::Index, const eic::ConstTrajectory*> rc2trajectory;
  for(const auto &trajectory: trajectories)
    rc2trajectory[trajectory.trackID()] = &trajectory;
#endif
    
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

    auto particle = new ChargedParticle(mctrack.pdgID());
    event->AddChargedParticle(particle);

    // Fill one global array of photon hits (shared by all radiators);
    std::vector<OpticalPhoton*> photons; 
    for(const auto &hit: hits) {
      // FIXME: yes, use MC truth here; not really needed I guess; 
      if (hit.g4ID() != mctrack.ID()) continue;
      
      // Simulate QE & geometric sensor efficiency; FIXME: hit.energy() is numerically 
      // in GeV units, but Gaudi::Units::GeV = 1000; prefer to convert photon energies 
      // to [eV] in all places by hand;
      if (!QE_pass(1E9*hit.energy(), m_rngUni()*m_GeometricEfficiency.value())) 
	continue;
      
      //printf("next photon() ... %5d %5d\n", hit.g4ID(), mctrack.ID());

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
      //unsigned module = (hit.cellID() >>  8) & 0x0FFF;
      photon->SetPhotonDetector(m_IrtDet->m_PhotonDetectors[0]);
      photon->SetDetected(true);
      photon->SetVolumeCopy(hit.cellID() & m_ReadoutCellMask);//module);
      
      photons.push_back(photon);
    } //for hit

    // Loop through the radiators one by one, using the same set of photons;
    for(auto radiator: m_SelectedRadiators) {
      // FIXME: yes, for now assume all the photons were produced in aerogel; 
      particle->StartRadiatorHistory(std::make_pair(radiator, new RadiatorHistory()));

      {
	// FIXME: hardcoded; give the algorithm radiator surface boundaries as encoded in ERich_geo.cpp;  
	auto s1 = radiator->GetFrontSide(0), s2 = radiator->GetRearSide(0);

	TVector3 from, to, p0;

#ifdef _USE_TRAJECTORIES_
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
	  const eic::TrajectoryPoint *best = 0;

	  // In case of aerogel, just find the closest point; gas case will be more complicated,
	  // but in general pretty much straightforward too;
	  for(const auto &point: trajectory->points()) 
	    // FIXME: optimize for efficiency later; now calculate every time new;
	    // Assume 's1' (front surface) is good enough;
	    if (!best || IRTAlgorithmServices::GetDistance(s1, &point) < 
		IRTAlgorithmServices::GetDistance(s1, best))
	      best = &point;

	  // FIXME: well, may want to check how far away the best point actually is;
	  IRTAlgorithmServices::GetCrossing(s1, best, &from);
	  IRTAlgorithmServices::GetCrossing(s2, best, &to);
	    
	  p0 = IRTAlgorithmServices::GetMomentum(best);
	} else {
#endif
	  // FIXME: need it not at vertex, but in the radiator; as coded here, this can 
	  // hardly work once the magnetic field is turned on; but this is a best guess 
	  // if nothing else is available;
	  const auto &vtx = mctrack.vs(), &p = mctrack.ps();
  
	  p0 = TVector3(p.x, p.y, p.z);
	  auto x0 = TVector3(vtx.x, vtx.y, vtx.z), n0 = p0.Unit();

	  s1->GetCrossing(x0, n0, &from);
	  s2->GetCrossing(x0, n0, &to);
	  
#ifdef _USE_TRAJECTORIES_
	} //if
#endif
	radiator->ResetLocations();

	{
	  unsigned zbins = radiator->GetTrajectoryBinCount();
	  TVector3 nn = (to - from).Unit(); from += (0.010)*nn; to -= (0.010)*nn;

	  // FIXME: assume a straight line to the moment;
	  auto span = to - from;
	  double tlen = span.Mag(), step = tlen / zbins;
	  for(unsigned iq=0; iq<zbins+1; iq++) 
	    radiator->AddLocation(from + iq*step*nn, p0);
	}
      }
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
      particle->PIDReconstruction(pid, &photons);
      
      // Now the interesting part comes: how to populate the output collection;
      auto aerogel = m_IrtDet->GetRadiator("Aerogel");
      for(auto radiator: m_SelectedRadiators) {
	// FIXME: extend the eicd table layout; FIXME: this loop is anyway inefficient;
	if (radiator != aerogel && m_SelectedRadiators.size() != 1) continue;
	
	auto cbuffer = cpid.create();
	
	for(unsigned ip=0; ip<sizeof(pdg_table)/sizeof(pdg_table[0]); ip++) {
	  auto hypo = pid.GetHypothesis(ip);
	  eic::CherenkovPdgHypothesis hypothesis;
	  
	  // Flip sign if needed; FIXME: use reconstructed one?;
	  if (m_pidSvc->particle(mctrack.pdgID()).charge < 0.0) pdg_table[ip] *= -1;
	  
	  hypothesis.pdg    = pdg_table[ip];
	  hypothesis.npe    = hypo->GetNpe   (radiator);//aerogel);
	  hypothesis.weight = hypo->GetWeight(radiator);//aerogel);
	  
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
      }
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


// Search a boolean solid's composition tree for primitive of type `typeName` (e.g., `TGeoSphere`)
// - `prim` will be set to the primitive; can be empty initially
// - `pos` will be set to the primitive's position (careful, it only considers translation of the operands);
//   should be initially set to the position of `sol`
// - `matx` stores operand transformations; you do not need to set it initially
void Jug::PID::IRTAlgorithm::findPrimitive(
    const std::string typeName, const dd4hep::Solid sol,
    dd4hep::Solid &prim, dd4hep::Position &pos, const TGeoMatrix *matx
) const {
  if(sol->IsComposite()) {
    auto node = (dd4hep::BooleanSolid) sol;
    findPrimitive(typeName, node.leftShape(),  prim, pos, node.leftMatrix()  );
    findPrimitive(typeName, node.rightShape(), prim, pos, node.rightMatrix() );
  } else if(typeName.compare(sol.type())==0) {
    prim = sol;
    //matx->Print();
    dd4hep::Position translation;
    translation.SetCoordinates(matx->GetTranslation());
    pos += translation;
  }
}

