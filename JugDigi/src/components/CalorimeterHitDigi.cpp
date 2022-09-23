// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Wouter Deconinck, Sylvester Joosten

#include <JugAlgo/Algorithm.h>
#include <algorithms/calorimetry/CalorimeterHitDigi.h>

#include "Gaudi/Property.h"

namespace {
  using AlgoBase = Jug::Algo::Algorithm<algorithms::calorimetry::CalorimeterHitDigi>;
}

class CalorimeterHitDigi : public AlgoBase {

public:
  CalorimeterHitDigi(const std::string& name, ISvcLocator* svcLoc) : AlgoBase(name, svcLoc) {}

  virtual StatusCode configure() {
    setAlgoProp("energyResolutions", u_eRes.value());
    setAlgoProp("timeResolution", m_tRes.value());
    setAlgoProp("threshold", m_threshold.value());
    setAlgoProp("capacityADC", m_capADC.value());
    setAlgoProp("dynamicRangeADC", m_dyRangeADC.value());
    setAlgoProp("pedMeanADC", m_pedMeanADC.value());
    setAlgoProp("pedSigmaADC", m_pedSigmaADC.value());
    setAlgoProp("resolutionTDC", m_resolutionTDC.value());
    setAlgoProp("scaleResponse", m_corrMeanScale.value());
    setAlgoProp("fieldRefNumbers", u_refs.value());
    setAlgoProp("geoServiceName", m_geoSvcName.value());
    setAlgoProp("readoutClass", m_readout.value());
    return StatusCode::SUCCESS;
  }

private:
  // additional smearing resolutions
  Gaudi::Property<std::vector<double>> u_eRes{this, "energyResolutions", {}}; // a/sqrt(E/GeV) + b + c/(E/GeV)
  Gaudi::Property<double>              m_tRes{this, "timeResolution", 0.0 * dd4hep::ns};
  // single hit energy deposition threshold
  Gaudi::Property<double>              m_threshold{this, "threshold", 1. * dd4hep::keV};

  // digitization settings
  Gaudi::Property<unsigned int>       m_capADC{this, "capacityADC", 8096};
  Gaudi::Property<double>             m_dyRangeADC{this, "dynamicRangeADC", 100 * dd4hep::MeV};
  Gaudi::Property<unsigned int>       m_pedMeanADC{this, "pedestalMean", 400};
  Gaudi::Property<double>             m_pedSigmaADC{this, "pedestalSigma", 3.2};
  Gaudi::Property<double>             m_resolutionTDC{this, "resolutionTDC", 0.010 * dd4hep::ns};
  Gaudi::Property<double>             m_corrMeanScale{this, "scaleResponse", 1.0};

  // signal sums
  // @TODO: implement signal sums with timing
  // field names to generate id mask, the hits will be grouped by masking the field
  Gaudi::Property<std::vector<std::string>> u_fields{this, "signalSumFields", {}};
  // ref field ids are used for the merged hits, 0 is used if nothing provided
  Gaudi::Property<std::vector<int>>         u_refs{this, "fieldRefNumbers", {}};
  Gaudi::Property<std::string>              m_geoSvcName{this, "geoServiceName", "GeoSvc"};
  Gaudi::Property<std::string>              m_readout{this, "readoutClass", ""};

};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(CalorimeterHitDigi)

} // namespace Jug::Digi
