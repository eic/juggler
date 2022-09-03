#pragma once

#include <functional>

// Algorithm base classes
#include "algorithms/base/algorithm.h"
#include "algorithms/base/property.h"
#include "algorithms/base/service.h"

// DD4hep
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDSegmentation/BitFieldCoder.h"

// EDM4hep
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"

namespace algorithms::digitization {

/** Generic calorimeter hit digitization.
 *
 * \ingroup digi
 * \ingroup calorimetry
 */
class CalorimeterHitDigi final : public JugAlgorithm<eicd::RawCalorimeterHitCollection, edm4hep::SimCalorimeterHitCollection> {

public:
  // additional smearing resolutions
  algorithms::Property<std::vector<double>> u_eRes{this, "energyResolutions", {}}; // a/sqrt(E/GeV) + b + c/(E/GeV)
  algorithms::Property<double> m_tRes{this, "timeResolution", 0.0 * dd4hep::ns};

  // digitization settings
  algorithms::Property<unsigned int> m_capADC{this, "capacityADC", 8096};
  algorithms::Property<double> m_dyRangeADC{this, "dynamicRangeADC", 100 * dd4hep::MeV};
  algorithms::Property<unsigned int> m_pedMeanADC{this, "pedestalMean", 400};
  algorithms::Property<double> m_pedSigmaADC{this, "pedestalSigma", 3.2};
  algorithms::Property<double> m_resolutionTDC{this, "resolutionTDC", 0.010 * dd4hep::ns};

  algorithms::Property<double> m_corrMeanScale{this, "scaleResponse", 1.0};
  // These are better variable names for the "energyResolutions" array which is a bit
  // magic @FIXME
  // algorithms::Property<double>             m_corrSigmaCoeffE{this, "responseCorrectionSigmaCoeffE", 0.0};
  // algorithms::Property<double>             m_corrSigmaCoeffSqrtE{this, "responseCorrectionSigmaCoeffSqrtE", 0.0};

  // signal sums
  // @TODO: implement signal sums with timing
  // field names to generate id mask, the hits will be grouped by masking the field
  algorithms::Property<std::vector<std::string>> u_fields{this, "signalSumFields", {}};
  // ref field ids are used for the merged hits, 0 is used if nothing provided
  algorithms::Property<std::vector<int>> u_refs{this, "fieldRefNumbers", {}};
  algorithms::Property<std::string> m_readout{this, "readoutClass", ""};

  // Geometry service
  algorithms::Service<dd4hep::Detector*(void)> m_geoSvc{this, "geoSvc"};

  // Random service
  algorithms::Service<double()> m_normDist{this, "normDist"};

  // unitless counterparts of inputs FIXME remove
  double dyRangeADC{0}, stepTDC{0}, tRes{0}, eRes[3] = {0., 0., 0.};
  uint64_t id_mask{0}, ref_mask{0};

  CalorimeterHitDigi() = default;

  bool initialize();

  eicd::RawCalorimeterHitCollection operator()(const edm4hep::SimCalorimeterHitCollection& input) const;

  bool finalize();

private:
  eicd::RawCalorimeterHitCollection single_hits_digi(const edm4hep::SimCalorimeterHitCollection& simhits) const;
  eicd::RawCalorimeterHitCollection signal_sum_digi(const edm4hep::SimCalorimeterHitCollection& simhits) const;
};

} // namespace algorithms::digitization
