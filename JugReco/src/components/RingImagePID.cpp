#include <algorithm>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/IParticlePropertySvc.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/TrackCollection.h"
#include "eicd/ReconstructedParticleCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  class RingImagePID : public GaudiAlgorithm {
  public:
    using Tracks       = eic::TrackCollection;
    using RecParticles = eic::ReconstructedParticleCollection;

    DataHandle<Tracks>     m_inputTracks{"inputTracks", Gaudi::DataHandle::Reader, this};
    DataHandle<RecParticles>    m_outputParticles{"outputParticles", Gaudi::DataHandle::Writer, this};
    Gaudi::Property<double> m_minModuleEdep{this, "minModuleEdep", 5.0 * MeV};
    Gaudi::Property<double> m_maxDistance{this, "maxDistance", 20.0 * cm};

    /// Pointer to the geometry service
    SmartIF<IGeoSvc>      m_geoSvc;
    IParticlePropertySvc* m_ppSvc;

    RingImagePID(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputTracks", m_inputTracks, "Input tracks without PID");
      declareProperty("outputParticles", m_outputParticles, "Output clusters");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      StatusCode sc = service("ParticlePropertySvc", m_ppSvc, true);
      if (sc.isFailure()) {
        error() << "Unable to locate Particle Property Service. "
                << "Check the configuration." << endmsg;
        return StatusCode::FAILURE;
      }
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override {
      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(RingImagePID)

} // namespace Jug::Reco
