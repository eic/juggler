// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Sylvester Joosten, Wouter Deconinck, Chao, Whitney Armstrong

// Reconstruct digitized outputs, paired with Jug::Digi::CalorimeterHitDigi
// Author: Chao Peng
// Date: 06/14/2021

#include <algorithms/algorithm.h>
#include <algorithms/geo.h>

// Event Model related classes
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/RawCalorimeterHitCollection.h"

namespace algorithms::calorimetry {

using CalorimeterHitRecoAlgorithm = Algorithm<
    Input<edm4eic::RawCalorimeterHitCollection>,
    Output<edm4eic::CalorimeterHitCollection>>;

/** Calorimeter hit reconstruction.
 *
 * Reconstruct digitized outputs, paired with Jug::Digi::CalorimeterHitDigi
 * \ingroup reco
 */
class CalorimeterHitReco : public CalorimeterHitRecoAlgorithm {

public:

  // TODO: get rid of "Collection" in names
  CalorimeterHitReco(std::string_view name)
      : CalorimeterHitRecoAlgorithm{name,
                            {"inputHitCollection"},
                            {"outputHitCollection"},
                            "Reconstruct digitized outputs, paired."} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  // length unit from dd4hep, should be fixed
  Property<double> m_lUnit{this, "lengthUnit", dd4hep::mm, "Length unit"};

  // digitization settings, must be consistent with digi class
  Property<unsigned int> m_capADC{this, "capacityADC", 8096, "Number of ADC channels"};
  Property<double> m_dyRangeADC{this, "dynamicRangeADC", 100. * dd4hep::MeV, "Dynamic range of the ADC"};
  Property<unsigned int> m_pedMeanADC{this, "pedestalMean", 400, "Mean of pedestal in ADC channels"};
  Property<double> m_pedSigmaADC{this, "pedestalSigma", 3.2, "Sigma of pedestal in ADC channels"};
  Property<double> m_resolutionTDC{this, "resolutionTDC", 10 * dd4hep::picosecond, ""};

  // zero suppression values
  Property<double> m_thresholdFactor{this, "thresholdFactor", 0.0, ""};
  Property<double> m_thresholdValue{this, "thresholdValue", 0.0, ""};

  // energy correction with sampling fraction
  Property<double> m_sampFrac{this, "samplingFraction", 1.0, ""};

  // unitless counterparts of the input parameters
  double dyRangeADC{0};
  double thresholdADC{0};
  double stepTDC{0};

  // geometry service to get ids, ignored if no names provided
  Property<std::string> m_geoSvcName{this, "geoServiceName", "GeoSvc"};
  Property<std::string> m_readout{this, "readoutClass", ""};
  Property<std::string> m_layerField{this, "layerField", ""};
  Property<std::string> m_sectorField{this, "sectorField", ""};

  dd4hep::BitFieldCoder* id_dec = nullptr;
  size_t sector_idx{0}, layer_idx{0};

  // name of detelment or fields to find the local detector (for global->local transform)
  // if nothing is provided, the lowest level DetElement (from cellID) will be used
  Property<std::string> m_localDetElement{this, "localDetElement", ""};
  Property<std::vector<std::string>> u_localDetFields{this, "localDetFields", {}};
  dd4hep::DetElement local;
  size_t local_mask = ~0;

  const GeoSvc& m_geo = GeoSvc::instance();
};
} // namespace algorithms::calorimetry
