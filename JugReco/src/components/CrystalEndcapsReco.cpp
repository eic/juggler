#include <algorithm>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {
  class CrystalEndcapsReco : public GaudiAlgorithm {
  public:
    Gaudi::Property<double>                      m_minModuleEdep{this, "minModuleEdep", 0.5 * MeV};
    DataHandle<eic::RawCalorimeterHitCollection> m_inputHitCollection{
        "inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection> m_outputHitCollection{
        "outputHitCollection", Gaudi::DataHandle::Writer, this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    CrystalEndcapsReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      warning() << "Deprecated algorithm for digi/reco, use Jug::Digi::CalorimeterHitDigi"
                   "and Jug::Reco::CalorimeterHitReco instead" << endmsg;
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& rawhits = *m_inputHitCollection.get();
      // Create output collections
      auto& hits = *m_outputHitCollection.createAndPut();

      // energy time reconstruction
      for (const auto& rh : rawhits) {
        float energy = rh.amplitude() / 1.0e6; // convert keV -> GeV
        if (energy >= (m_minModuleEdep / GeV)) {
          float time = rh.timeStamp();
          auto  id   = rh.cellID();
          // global positions
          auto gpos = m_geoSvc->cellIDPositionConverter()->position(id);
          // local positions
          auto pos =
              m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().position();
          // cell dimension
          auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(id);
          hits.push_back(eic::CalorimeterHit{id,
                                             -1,
                                             -1,
                                             -1,
                                             energy,
                                             time,
                                             {gpos.x(), gpos.y(), gpos.z()},
                                             {pos.x(), pos.y(), pos.z()},
                                             {dim[0], dim[1], 0.0},
                                             0});
        }
      }

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(CrystalEndcapsReco)

} // namespace Jug::Reco
