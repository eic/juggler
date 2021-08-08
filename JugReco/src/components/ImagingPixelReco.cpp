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

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"

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
  public:
    // geometry service
    Gaudi::Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
    Gaudi::Property<std::string> m_readout{this, "readoutClass", ""};
    Gaudi::Property<std::string> m_layerField{this, "layerField", "layer"};
    Gaudi::Property<std::string> m_sectorField{this, "sectorField", "sector"};
    // length unit (from dd4hep geometry service)
    Gaudi::Property<double> m_lUnit{this, "lengthUnit", dd4hep::mm};
    // digitization parameters
    Gaudi::Property<int>    m_capADC{this, "capacityADC", 8096};
    Gaudi::Property<int>    m_pedMeanADC{this, "pedestalMean", 400};
    Gaudi::Property<double> m_dyRangeADC{this, "dynamicRangeADC", 100 * MeV};
    Gaudi::Property<double> m_pedSigmaADC{this, "pedestalSigma", 3.2};
    Gaudi::Property<double> m_thresholdADC{this, "thresholdFactor", 3.0};
    // Calibration!
    Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};

    // unitless counterparts for the input parameters
    double dyRangeADC;

    // hits containers
    DataHandle<eic::RawCalorimeterHitCollection> m_inputHitCollection{
        "inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection> m_outputHitCollection{"outputHitCollection",
                                                                    Gaudi::DataHandle::Writer, this};

    // Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;
    // visit readout fields
    dd4hep::BitFieldCoder* id_dec;
    size_t                 sector_idx, layer_idx;

    ImagingPixelReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
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

      if (m_readout.value().empty()) {
        error() << "readoutClass is not provided, it is needed to know the fields in readout ids"
                << endmsg;
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
      dyRangeADC = m_dyRangeADC.value()/GeV;

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& rawhits = *m_inputHitCollection.get();
      // Create output collections
      auto& hits = *m_outputHitCollection.createAndPut();

      // energy time reconstruction
      int nhits = 0;
      for (const auto& rh : rawhits) {
        // did not pass the threshold
        if ((rh.amplitude() - m_pedMeanADC) < m_thresholdADC * m_pedSigmaADC) {
          continue;
        }
        double energy = (rh.amplitude() - m_pedMeanADC) / (double)m_capADC * dyRangeADC / m_sampFrac;   // convert ADC -> energy
        double time = rh.time() * 1.e-6; // ns
        auto   id   = rh.cellID();
        int    lid  = (int)id_dec->get(id, layer_idx);
        int    sid  = (int)id_dec->get(id, sector_idx);

        // global positions
        auto gpos = m_geoSvc->cellIDPositionConverter()->position(id);
        // local positions
        auto volman    = m_geoSvc->detector()->volumeManager();
        auto alignment = volman.lookupDetElement(id).nominal();
        auto pos       = alignment.worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));

        hits.push_back(eic::CalorimeterHit{
            id,           // cellID
            nhits++,      // ID
            lid,          // layer
            sid,          // sector
            0,            // type
            static_cast<float>(energy), // energy
            0,                          // energyError
            static_cast<float>(time),   // time
            {gpos.x() / m_lUnit, gpos.y() / m_lUnit, gpos.z() / m_lUnit}, // global pos
            {pos.x() / m_lUnit, pos.y() / m_lUnit, pos.z() / m_lUnit},    // local pos
            {0, 0, 0}}); // @TODO: add dimension
      }
      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(ImagingPixelReco)

} // namespace Jug::Reco
