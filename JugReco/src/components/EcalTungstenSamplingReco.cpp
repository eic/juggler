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
    Gaudi::Property<double>                      m_lUnit{this, "lengthUnit", dd4hep::mm};
    Gaudi::Property<int>                         m_capADC{this, "capacityADC", 8096};
    Gaudi::Property<double>                      m_dyRangeADC{this, "dynamicRangeADC", 100 * MeV};
    Gaudi::Property<int>                         m_pedMeanADC{this, "pedestalMean", 400};
    Gaudi::Property<double>                      m_pedSigmaADC{this, "pedestalSigma", 3.2};
    Gaudi::Property<double>                      m_thresholdADC{this, "thresholdFactor", 3.0};
    DataHandle<eic::RawCalorimeterHitCollection> m_inputHitCollection{
        "inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection> m_outputHitCollection{
        "outputHitCollection", Gaudi::DataHandle::Writer, this};
    // geometry service
    Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
    Gaudi::Property<std::string> m_readout{this, "readoutClass", ""};
    Gaudi::Property<std::string> m_layerField{this, "layerField", ""};
    Gaudi::Property<std::string> m_sectorField{this, "sectorField", ""};
    SmartIF<IGeoSvc>             m_geoSvc;
    dd4hep::BitFieldCoder*       id_dec = nullptr;
    size_t                       sector_idx, layer_idx;

    EcalTungstenSamplingReco(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
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

      m_geoSvc = service(m_geoSvcName);
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }

      // do not get the layer/sector ID if no readout class provided
      if (m_readout.value().empty()) {
        return StatusCode::SUCCESS;
      }

      try {
        id_dec = m_geoSvc->detector()->readout(m_readout).idSpec().decoder();
        if (m_sectorField.value().size()) {
          sector_idx = id_dec->index(m_sectorField);
        }
        if (m_layerField.value().size()) {
          layer_idx = id_dec->index(m_layerField);
        }
      } catch (...) {
        error() << "Failed to load ID decoder for " << m_readout << endmsg;
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
        // did not pass the threshold
        if ((rh.amplitude() - m_pedMeanADC) < m_thresholdADC * m_pedSigmaADC) {
          continue;
        }
        float energy = (rh.amplitude() - m_pedMeanADC) / (float)m_capADC *
                       m_dyRangeADC; // convert ADC -> energy
        float time = rh.timeStamp(); // ns
        auto  id   = rh.cellID();
        int   lid  = ((id_dec != nullptr) & m_layerField.value().size())
                      ? static_cast<int>(id_dec->get(id, layer_idx))
                      : -1;
        int sid = ((id_dec != nullptr) & m_sectorField.value().size())
                      ? static_cast<int>(id_dec->get(id, sector_idx))
                      : -1;

        // global positions
        auto gpos = m_geoSvc->cellIDPositionConverter()->position(id);
        // local positions
        auto volman    = m_geoSvc->detector()->volumeManager();
        auto alignment = volman.lookupDetElement(id).nominal();
        auto pos       = alignment.worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
        // auto pos =
        // m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().position(); cell
        // dimension
        auto   cdim   = m_geoSvc->cellIDPositionConverter()->cellDimensions(id);
        double dim[3] = {0., 0., 0.};
        for (size_t i = 0; i < cdim.size() && i < 3; ++i) {
          dim[i] = cdim[i] / m_lUnit;
        }
        // info() << std::bitset<64>(id) << "\n"
        //        <<
        //        m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().volIDs().str()
        //        << endmsg;
        hits.push_back(
            eic::CalorimeterHit{id,
                                -1,
                                lid,
                                sid,
                                energy,
                                time,
                                {gpos.x() / m_lUnit, gpos.y() / m_lUnit, gpos.z() / m_lUnit},
                                {pos.x() / m_lUnit, pos.y() / m_lUnit, pos.z() / m_lUnit},
                                {dim[0], dim[1], dim[2]},
                                0});
      }

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(EcalTungstenSamplingReco)

} // namespace Jug::Reco
