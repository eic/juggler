// Reconstruct digitized outputs, paired with Jug::Digi::CalorimeterHitDigi
// Author: Chao Peng
// Date: 06/14/2021

#include "fmt/format.h"
#include "fmt/ranges.h"
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

  class CalorimeterHitReco : public GaudiAlgorithm {
  public:
    // length unit from dd4hep, should be fixed
    Gaudi::Property<double> m_lUnit{this, "lengthUnit", dd4hep::mm};

    // digitization settings, must be consistent with digi class
    Gaudi::Property<int>    m_capADC{this, "capacityADC", 8096};
    Gaudi::Property<double> m_dyRangeADC{this, "dynamicRangeADC", 100. * MeV};
    Gaudi::Property<int>    m_pedMeanADC{this, "pedestalMean", 400};
    Gaudi::Property<double> m_pedSigmaADC{this, "pedestalSigma", 3.2};

    // zero suppression values
    Gaudi::Property<double> m_thresholdADC{this, "thresholdFactor", 3.0};

    // unitless counterparts of the input parameters
    double                  dyRangeADC;

    DataHandle<eic::RawCalorimeterHitCollection> m_inputHitCollection{
        "inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection> m_outputHitCollection{
        "outputHitCollection", Gaudi::DataHandle::Writer, this};

    // geometry service to get ids, ignored if no names provided
    Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
    Gaudi::Property<std::string> m_readout{this, "readoutClass", ""};
    Gaudi::Property<std::string> m_layerField{this, "layerField", ""};
    Gaudi::Property<std::string> m_sectorField{this, "sectorField", ""};
    SmartIF<IGeoSvc>             m_geoSvc;
    dd4hep::BitFieldCoder*       id_dec = nullptr;
    size_t                       sector_idx, layer_idx;

    // name of detelment or fields to find the local detector (for global->local transform)
    // if nothing is provided, the lowest level DetElement (from cellID) will be used
    Gaudi::Property<std::string>              m_localDetElement{this, "localDetElement", ""};
    Gaudi::Property<std::vector<std::string>> u_localDetFields{this, "localDetFields", {}};
    dd4hep::DetElement                        local;
    size_t                                    local_mask = ~0;

    CalorimeterHitReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
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

      // unitless conversion
      dyRangeADC = m_dyRangeADC.value()/GeV;

      // do not get the layer/sector ID if no readout class provided
      if (m_readout.value().empty()) {
        return StatusCode::SUCCESS;
      }

      auto id_spec = m_geoSvc->detector()->readout(m_readout).idSpec();
      try {
        id_dec = id_spec.decoder();
        if (m_sectorField.value().size()) {
          sector_idx = id_dec->index(m_sectorField);
          info() << "Find sector field " << m_sectorField.value() << ", index = " << sector_idx << endmsg;
        }
        if (m_layerField.value().size()) {
          layer_idx = id_dec->index(m_layerField);
          info() << "Find layer field " << m_layerField.value() << ", index = " << sector_idx << endmsg;
        }
      } catch (...) {
        error() << "Failed to load ID decoder for " << m_readout << endmsg;
        return StatusCode::FAILURE;
      }

      // local detector name has higher priority
      if (m_localDetElement.value().size()) {
        try {
          local = m_geoSvc->detector()->detector(m_localDetElement.value());
          info() << "Local coordinate system from DetElement " << m_localDetElement.value()
                 << endmsg;
        } catch (...) {
          error() << "Failed to locate local coordinate system from DetElement "
                  << m_localDetElement.value() << endmsg;
          return StatusCode::FAILURE;
        }
      // or get from fields
      } else {
        std::vector<std::pair<std::string, int>> fields;
        for (auto& f : u_localDetFields.value()) {
          fields.push_back({f, 0});
        }
        local_mask = id_spec.get_mask(fields);
        // use all fields if nothing provided
        if (fields.empty()) {
          local_mask = ~0;
        }
        info() << fmt::format("Local DetElement mask {:#064b} from fields [{}]", local_mask,
                              fmt::join(fields, ", "))
               << endmsg;
      }

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& rawhits = *m_inputHitCollection.get();
      // create output collections
      auto& hits = *m_outputHitCollection.createAndPut();
      auto converter = m_geoSvc->cellIDPositionConverter();

      // energy time reconstruction
      for (const auto& rh : rawhits) {
        // did not pass the zero-suppression threshold
        if ((rh.amplitude() - m_pedMeanADC) < m_thresholdADC * m_pedSigmaADC) {
          continue;
        }

        // convert ADC -> energy
        float energy = (rh.amplitude() - m_pedMeanADC) / static_cast<float>(m_capADC.value()) * dyRangeADC;

        float time = rh.timeStamp(); // ns
        auto  id   = rh.cellID();
        int   lid  = ((id_dec != nullptr) && m_layerField.value().size())
                      ? static_cast<int>(id_dec->get(id, layer_idx))
                      : -1;
        int   sid  = ((id_dec != nullptr) && m_sectorField.value().size())
                      ? static_cast<int>(id_dec->get(id, sector_idx))
                      : -1;
        // global positions
        auto gpos = converter->position(id);
        // local positions
        if (m_localDetElement.value().empty()) {
          auto volman = m_geoSvc->detector()->volumeManager();
          local       = volman.lookupDetElement(id & local_mask);
        }
        auto pos = local.nominal().worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
        // auto pos = m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().position();
        // cell dimension
        std::vector<double> cdim;
        // get segmentation dimensions
        if (converter->findReadout(local).segmentation().type() != "NoSegmentation") {
          cdim  = converter->cellDimensions(id);
        // get volume dimensions (multiply by two to get fullsize)
        } else {
          // cdim = converter->findContext(id)->volumePlacement().volume().solid().dimensions();
          // Using bounding box instead of actual solid so the dimensions are always in dim_x, dim_y, dim_z
          cdim = converter->findContext(id)->volumePlacement().volume().boundingBox().dimensions();
          std::transform(cdim.begin(), cdim.end(), cdim.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, 2));
        }
        double dim[3] = {0., 0., 0.};
        for (size_t i = 0; i < cdim.size() && i < 3; ++i) {
          dim[i] = cdim[i] / m_lUnit;
        }
        // debug() << std::bitset<64>(id) << "\n"
        //         <<
        //         m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().volIDs().str()
        //         << endmsg;
        hits.push_back({
            id,
            -1,
            lid,
            sid, // cell id, cluster id, layer id, sector id
            energy,
            time, // energy, time
            {gpos.x() / m_lUnit, gpos.y() / m_lUnit, gpos.z() / m_lUnit},
            {pos.x() / m_lUnit, pos.y() / m_lUnit, pos.z() / m_lUnit},
            {dim[0], dim[1], dim[2]},
            0 // @TODO: hit type
        });
      }

      return StatusCode::SUCCESS;
    }

  }; // class CalorimeterHitReco

  DECLARE_COMPONENT(CalorimeterHitReco)

} // namespace Jug::Reco
