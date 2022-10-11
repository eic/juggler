// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Sylvester Joosten, Chao, Whitney Armstrong, Wouter Deconinck, Jihee Kim

// Reconstruct digitized outputs of ImagingCalorimeter
// It converts digitized ADC/TDC values to energy/time, and looks for geometrical information of the
// readout pixels Author: Chao Peng Date: 06/02/2021

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
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/RawCalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

/** Imaging calorimeter pixel hit reconstruction.
 *
 * Reconstruct digitized outputs of ImagingCalorimeter
 * It converts digitized ADC/TDC values to energy/time, and looks for geometrical information of the
 *
 * \ingroup reco
 */
class ImagingPixelReco : public GaudiAlgorithm {
private:
  // geometry service
  Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
  Gaudi::Property<std::string> m_readout{this, "readoutClass", ""};
  Gaudi::Property<std::string> m_layerField{this, "layerField", "layer"};
  Gaudi::Property<std::string> m_sectorField{this, "sectorField", "sector"};
  // length unit (from dd4hep geometry service)
  Gaudi::Property<double> m_lUnit{this, "lengthUnit", dd4hep::mm};
  // digitization parameters
  Gaudi::Property<unsigned int> m_capADC{this, "capacityADC", 8096};
  Gaudi::Property<unsigned int> m_pedMeanADC{this, "pedestalMean", 400};
  Gaudi::Property<double> m_dyRangeADC{this, "dynamicRangeADC", 100 * MeV};
  Gaudi::Property<double> m_pedSigmaADC{this, "pedestalSigma", 3.2};
  Gaudi::Property<double> m_thresholdADC{this, "thresholdFactor", 3.0};
  // Calibration!
  Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};

  // unitless counterparts for the input parameters
  double dyRangeADC{0};

  // hits containers
  DataHandle<edm4eic::RawCalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                    this};
  DataHandle<edm4eic::CalorimeterHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer,
                                                                  this};

  // Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  // visit readout fields
  dd4hep::BitFieldCoder* id_dec;
  size_t sector_idx{0}, layer_idx{0};

public:
  ImagingPixelReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
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

    if (m_readout.value().empty()) {
      error() << "readoutClass is not provided, it is needed to know the fields in readout ids" << endmsg;
      return StatusCode::FAILURE;
    }

    try {
      id_dec     = m_geoSvc->detector()->readout(m_readout).idSpec().decoder();
      sector_idx = id_dec->index(m_sectorField);
      layer_idx  = id_dec->index(m_layerField);
    } catch (...) {
      error() << "Failed to load ID decoder for " << m_readout << endmsg;
      return StatusCode::FAILURE;
    }

    // unitless conversion
    dyRangeADC = m_dyRangeADC.value() / GeV;

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& rawhits = *m_inputHitCollection.get();
    // Create output collections
    auto& hits = *m_outputHitCollection.createAndPut();

    // energy time reconstruction
    for (const auto& rh : rawhits) {

      #pragma GCC diagnostic push
      #pragma GCC diagnostic error "-Wsign-conversion"

      // did not pass the threshold
      if (rh.getAmplitude() < m_pedMeanADC + m_thresholdADC * m_pedSigmaADC) {
        continue;
      }
      const double energy =
        (((signed)rh.getAmplitude() - (signed)m_pedMeanADC)) / (double)m_capADC * dyRangeADC / m_sampFrac; // convert ADC -> energy
      const double time = rh.getTimeStamp() * 1.e-6;                                       // ns

      #pragma GCC diagnostic pop

      const auto id     = rh.getCellID();
      // @TODO remove
      const int lid = (int)id_dec->get(id, layer_idx);
      const int sid = (int)id_dec->get(id, sector_idx);

      // global positions
      const auto gpos = m_geoSvc->cellIDPositionConverter()->position(id);
      // local positions
      const auto volman = m_geoSvc->detector()->volumeManager();
      // TODO remove
      const auto alignment = volman.lookupDetElement(id).nominal();
      const auto pos       = alignment.worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));


      // create const vectors for passing to hit initializer list
      const decltype(edm4eic::CalorimeterHitData::position) position(
        gpos.x() / m_lUnit, gpos.y() / m_lUnit, gpos.z() / m_lUnit
      );
      const decltype(edm4eic::CalorimeterHitData::local) local(
        pos.x() / m_lUnit, pos.y() / m_lUnit, pos.z() / m_lUnit
      );

      hits.push_back(edm4eic::CalorimeterHit{id,                         // cellID
                                          static_cast<float>(energy), // energy
                                          0,                          // energyError
                                          static_cast<float>(time),   // time
                                          0,                          // timeError TODO
                                          position,                   // global pos
                                          {0, 0, 0}, // @TODO: add dimension
                                          sid,lid,
                                          local});                    // local pos
    }
    return StatusCode::SUCCESS;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ImagingPixelReco)

} // namespace Jug::Reco
