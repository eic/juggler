// Reconstruct digitized outputs fof Ecal Tungsten Sampling Calorimeter
// It is exactly the reverse step of JugDigi/src/components/EcalTungstenSamplingDigi.cpp

#include <algorithm>
#include <bitset>

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
  class EcalTungstenSamplingReco : public GaudiAlgorithm {
  public:
    Gaudi::Property<int>                         m_capADC{this, "capacityADC", 8096};
    Gaudi::Property<double>                      m_dyRangeADC{this, "DynamicRangeADC", 100*MeV};
    Gaudi::Property<int>                         m_pedMeanADC{this, "pedestalMean", 400};
    Gaudi::Property<double>                      m_pedSigmaADC{this, "pedestalSigma", 3.2};
    Gaudi::Property<double>                      m_thresholdADC{this, "thresholdFactor", 3.0};
    DataHandle<eic::RawCalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                      this};
    DataHandle<eic::CalorimeterHitCollection>    m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                    this};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    EcalTungstenSamplingReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
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

    StatusCode execute() override
    {
      // input collections
      const auto& rawhits = *m_inputHitCollection.get();
      // Create output collections
      auto& hits = *m_outputHitCollection.createAndPut();

      // energy time reconstruction
      for (auto& rh : rawhits) {
        // did not pass the threshold
        if ((rh.amplitude() - m_pedMeanADC) < m_thresholdADC*m_pedSigmaADC) {
          continue;
        }
        float energy = (rh.amplitude() - m_pedMeanADC) / (float) m_capADC * m_dyRangeADC; // convert ADC -> energy
        float time = rh.timeStamp(); // ns
        auto id = rh.cellID();
        // global positions
        auto gpos = m_geoSvc->cellIDPositionConverter()->position(id);
        // local positions
        auto volman = m_geoSvc->detector()->volumeManager();
        auto alignment = volman.lookupDetector(id).nominal();
        auto pos = alignment.worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
        // auto pos = m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().position();
        // cell dimension
        auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(id);
        // info() << std::bitset<64>(id) << "\n"
        //        << m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().volIDs().str() << endmsg;
        hits.push_back(eic::CalorimeterHit{id,
                                           energy,
                                           time,
                                           {gpos.x() / dd4hep::mm, gpos.y() / dd4hep::mm, gpos.z() / dd4hep::mm},
                                           {pos.x() / dd4hep::mm, pos.y() / dd4hep::mm, pos.z() / dd4hep::mm},
                                           {dim[0] / dd4hep::mm, dim[1] / dd4hep::mm, 0.0},
                                           0});
      }

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(EcalTungstenSamplingReco)

} // namespace Jug::Reco
