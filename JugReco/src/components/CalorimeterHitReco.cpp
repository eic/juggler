// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Sylvester Joosten, Wouter Deconinck, Chao, Whitney Armstrong

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

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

/** Calorimeter hit reconstruction.
 *
 * Reconstruct digitized outputs, paired with Jug::Digi::CalorimeterHitDigi
 * \ingroup reco
 */
class CalorimeterHitReco : public GaudiAlgorithm {
private:
  // length unit from dd4hep, should be fixed
  Gaudi::Property<double> m_lUnit{this, "lengthUnit", dd4hep::mm};

  // digitization settings, must be consistent with digi class
  Gaudi::Property<int> m_capADC{this, "capacityADC", 8096};
  Gaudi::Property<double> m_dyRangeADC{this, "dynamicRangeADC", 100. * MeV};
  Gaudi::Property<int> m_pedMeanADC{this, "pedestalMean", 400};
  Gaudi::Property<double> m_pedSigmaADC{this, "pedestalSigma", 3.2};
  Gaudi::Property<double> m_resolutionTDC{this, "resolutionTDC", 10 * ps};

  // zero suppression values
  Gaudi::Property<double> m_thresholdFactor{this, "thresholdFactor", 0.0};
  Gaudi::Property<double> m_thresholdValue{this, "thresholdValue", 0.0};

  // energy correction with sampling fraction
  Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};

  // unitless counterparts of the input parameters
  double dyRangeADC{0};
  double thresholdADC{0};
  double stepTDC{0};

  DataHandle<eicd::RawCalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                    this};
  DataHandle<eicd::CalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                  this};

  // geometry service to get ids, ignored if no names provided
  Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
  Gaudi::Property<std::string> m_readout{this, "readoutClass", ""};
  Gaudi::Property<std::string> m_layerField{this, "layerField", ""};
  Gaudi::Property<std::string> m_sectorField{this, "sectorField", ""};
  SmartIF<IGeoSvc> m_geoSvc;
  dd4hep::BitFieldCoder* id_dec = nullptr;
  size_t sector_idx{0}, layer_idx{0};

  // name of detelment or fields to find the local detector (for global->local transform)
  // if nothing is provided, the lowest level DetElement (from cellID) will be used
  Gaudi::Property<std::string> m_localDetElement{this, "localDetElement", ""};
  Gaudi::Property<std::vector<std::string>> u_localDetFields{this, "localDetFields", {}};
  dd4hep::DetElement local;
  size_t local_mask = ~0;

public:
  CalorimeterHitReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputHitCollection", m_outputHitCollection, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    m_geoSvc = service(m_geoSvcName);
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }

    // unitless conversion
    dyRangeADC = m_dyRangeADC.value() / GeV;

    // threshold for firing
    thresholdADC = m_thresholdFactor.value() * m_pedSigmaADC.value() + m_thresholdValue.value();

    // TDC channels to timing conversion
    stepTDC = ns / m_resolutionTDC.value();

    // do not get the layer/sector ID if no readout class provided
    if (m_readout.value().empty()) {
      return StatusCode::SUCCESS;
    }

    auto id_spec = m_geoSvc->detector()->readout(m_readout).idSpec();
    try {
      id_dec = id_spec.decoder();
      if (!m_sectorField.value().empty()) {
        sector_idx = id_dec->index(m_sectorField);
        info() << "Find sector field " << m_sectorField.value() << ", index = " << sector_idx << endmsg;
      }
      if (!m_layerField.value().empty()) {
        layer_idx = id_dec->index(m_layerField);
        info() << "Find layer field " << m_layerField.value() << ", index = " << sector_idx << endmsg;
      }
    } catch (...) {
      error() << "Failed to load ID decoder for " << m_readout << endmsg;
      return StatusCode::FAILURE;
    }

    // local detector name has higher priority
    if (!m_localDetElement.value().empty()) {
      try {
        local = m_geoSvc->detector()->detector(m_localDetElement.value());
        info() << "Local coordinate system from DetElement " << m_localDetElement.value() << endmsg;
      } catch (...) {
        error() << "Failed to locate local coordinate system from DetElement " << m_localDetElement.value() << endmsg;
        return StatusCode::FAILURE;
      }
      // or get from fields
    } else {
      std::vector<std::pair<std::string, int>> fields;
      for (auto& f : u_localDetFields.value()) {
        fields.emplace_back(f, 0);
      }
      local_mask = id_spec.get_mask(fields);
      // use all fields if nothing provided
      if (fields.empty()) {
        local_mask = ~0;
      }
      info() << fmt::format("Local DetElement mask {:#064b} from fields [{}]", local_mask, fmt::join(fields, ", "))
             << endmsg;
    }

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& rawhits = *m_inputHitCollection.get();
    // create output collections
    auto& hits     = *m_outputHitCollection.createAndPut();
    auto converter = m_geoSvc->cellIDPositionConverter();

    // energy time reconstruction
    for (const auto& rh : rawhits) {
      // did not pass the zero-suppression threshold
      if ((signed)(rh.getAmplitude() - m_pedMeanADC) < thresholdADC) {
        continue;
      }

      // convert ADC -> energy
      const float energy =
          (signed)(rh.getAmplitude() - m_pedMeanADC) / static_cast<float>(m_capADC.value()) * dyRangeADC / m_sampFrac; 
      const float time  = rh.getTimeStamp() / stepTDC;
      const auto cellID = rh.getCellID();
      const int lid = id_dec != nullptr && !m_layerField.value().empty() ? static_cast<int>(id_dec->get(cellID, layer_idx)) : -1;
      const int sid = id_dec != nullptr && !m_sectorField.value().empty() ? static_cast<int>(id_dec->get(cellID, sector_idx)) : -1;
      // global positions
      const auto gpos = converter->position(cellID);
      // local positions
      if (m_localDetElement.value().empty()) {
        auto volman = m_geoSvc->detector()->volumeManager();
        local       = volman.lookupDetElement(cellID & local_mask);
        }
        const auto pos = local.nominal().worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
        // cell dimension
        std::vector<double> cdim;
        // get segmentation dimensions
        if (converter->findReadout(local).segmentation().type() != "NoSegmentation") {
        cdim = converter->cellDimensions(cellID);
        // get volume dimensions (multiply by two to get fullsize)
        } else {
        // Using bounding box instead of actual solid so the dimensions are always in dim_x, dim_y, dim_z
        cdim = converter->findContext(cellID)->volumePlacement().volume().boundingBox().dimensions();
        std::transform(cdim.begin(), cdim.end(), cdim.begin(),
                       std::bind(std::multiplies<double>(), std::placeholders::_1, 2));
        }

        // create const vectors for passing to hit initializer list
        const decltype(eicd::CalorimeterHitData::position) position(
          gpos.x() / m_lUnit, gpos.y() / m_lUnit, gpos.z() / m_lUnit
        );
        const decltype(eicd::CalorimeterHitData::dimension) dimension(
          cdim[0] / m_lUnit, cdim[1] / m_lUnit, cdim[2] / m_lUnit
        );
        const decltype(eicd::CalorimeterHitData::local) local_position(
          pos.x() / m_lUnit, pos.y() / m_lUnit, pos.z() / m_lUnit
        );

        hits.push_back({
            rh.getCellID(), // cellID
            energy,         // energy
            0,              // @TODO: energy error
            time,           // time
            0,              // time error FIXME should be configurable
            position,       // global pos
            dimension,
            // Local hit info
            sid,
            lid,
            local_position, // local pos
        });
    }

    return StatusCode::SUCCESS;
  }

}; // class CalorimeterHitReco

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(CalorimeterHitReco)

} // namespace Jug::Reco
