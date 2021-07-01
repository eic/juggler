#include <algorithm>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DD4hep/DD4hepUnits.h"
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

  class EMCalReconstruction : public GaudiAlgorithm {
  public:
    using RawHits = eic::RawCalorimeterHitCollection;
    using RecHits = eic::CalorimeterHitCollection;

    DataHandle<RawHits> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<RecHits> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                              this};
    Gaudi::Property<double> m_minModuleEdep{this, "minModuleEdep", 5.0 * MeV};
    Gaudi::Property<double> m_samplingFraction{this, "samplingFraction", 1.0};

    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    EMCalReconstruction(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      warning() << "Deprecated algorithm for digi/reco, use Jug::Digi::CalorimeterHitDigi"
                   "and Jug::Reco::CalorimeterHitReco instead"
                << endmsg;
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
        float energy = rh.amplitude() / 1e6; // GeV
        if (energy >= m_minModuleEdep) {
          float time = rh.timeStamp() / 1e6; // ns
          auto  id   = rh.cellID();
          // global positions
          auto gpos = m_geoSvc->cellIDPositionConverter()->position(id);
          // local positions
          auto pos =
              m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().position();
          // cell dimension
          auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(id);
          hits.push_back(eic::CalorimeterHit{
              id,
              -1,
              -1,
              -1,
              static_cast<float>(energy / m_samplingFraction),
              time,
              {gpos.x() / dd4hep::mm, gpos.y() / dd4hep::mm, gpos.z() / dd4hep::mm},
              {pos.x() / dd4hep::mm, pos.y() / dd4hep::mm, pos.z() / dd4hep::mm},
              {dim[0] / dd4hep::mm, dim[1] / dd4hep::mm, 0.0},
              0});
        }
      }

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(EMCalReconstruction)

} // namespace Jug::Reco
