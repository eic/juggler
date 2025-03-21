// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Wouter Deconinck, Sylvester Joosten

// A general digitization for CalorimeterHit from simulation
// 1. Smear energy deposit with a/sqrt(E/GeV) + b + c/E or a/sqrt(E/GeV) (relative value)
// 2. Digitize the energy with dynamic ADC range and add pedestal (mean +- sigma)
// 3. Time conversion with smearing resolution (absolute value)
// 4. Signal is summed if the SumFields are provided
//
// Author: Chao Peng
// Date: 06/02/2021

#include <algorithm>
#include <cmath>
#include <unordered_map>

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/RndmGenerators.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDSegmentation/BitFieldCoder.h"

#include <k4Interface/IGeoSvc.h>
#include <k4FWCore/DataHandle.h>

#include "fmt/format.h"
#include "fmt/ranges.h"

// Event Model related classes
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"


using namespace Gaudi::Units;

namespace Jug::Digi {

  /** Generic calorimeter hit digitiziation.
   *
   * \ingroup digi
   * \ingroup calorimetry
   */
  class CalorimeterHitDigi : public Gaudi::Algorithm {
  public:
    // additional smearing resolutions
    Gaudi::Property<std::vector<double>> u_eRes{this, "energyResolutions", {}}; // a/sqrt(E/GeV) + b + c/(E/GeV)
    Gaudi::Property<double>              m_tRes{this, "timeResolution", 0.0 * ns};
    // single hit energy deposition threshold
    Gaudi::Property<double>              m_threshold{this, "threshold", 1. * keV};

    // digitization settings
    Gaudi::Property<unsigned int>       m_capADC{this, "capacityADC", 8096};
    Gaudi::Property<double>             m_dyRangeADC{this, "dynamicRangeADC", 100 * MeV};
    Gaudi::Property<unsigned int>       m_pedMeanADC{this, "pedestalMean", 400};
    Gaudi::Property<double>             m_pedSigmaADC{this, "pedestalSigma", 3.2};
    Gaudi::Property<double>             m_resolutionTDC{this, "resolutionTDC", 10 * ps};

    Gaudi::Property<double>             m_corrMeanScale{this, "scaleResponse", 1.0};
    // These are better variable names for the "energyResolutions" array which is a bit
    // magic @FIXME
    //Gaudi::Property<double>             m_corrSigmaCoeffE{this, "responseCorrectionSigmaCoeffE", 0.0};
    //Gaudi::Property<double>             m_corrSigmaCoeffSqrtE{this, "responseCorrectionSigmaCoeffSqrtE", 0.0};

    // signal sums
    // @TODO: implement signal sums with timing
    // field names to generate id mask, the hits will be grouped by masking the field
    Gaudi::Property<std::vector<std::string>> u_fields{this, "signalSumFields", {}};
    // ref field ids are used for the merged hits, 0 is used if nothing provided
    Gaudi::Property<std::vector<int>>         u_refs{this, "fieldRefNumbers", {}};
    Gaudi::Property<std::string>              m_geoSvcName{this, "geoServiceName", "GeoSvc"};
    Gaudi::Property<std::string>              m_readout{this, "readoutClass", ""};

    // unitless counterparts of inputs
    double           dyRangeADC{0}, stepTDC{0}, tRes{0}, eRes[3] = {0., 0., 0.};
    Rndm::Numbers    m_normDist;
    SmartIF<IGeoSvc> m_geoSvc;
    uint64_t         id_mask{0}, ref_mask{0};

    mutable DataHandle<edm4hep::SimCalorimeterHitCollection> m_inputHitCollection{
      "inputHitCollection", Gaudi::DataHandle::Reader, this};
    mutable DataHandle<edm4hep::RawCalorimeterHitCollection> m_outputHitCollection{
      "outputHitCollection", Gaudi::DataHandle::Writer, this};

    //  ill-formed: using Gaudi::Algorithm::GaudiAlgorithm;
    CalorimeterHitDigi(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      if (Gaudi::Algorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      // random number generator from service
      auto randSvc = Gaudi::svcLocator()->service<IRndmGenSvc>("RndmGenSvc", true);
      auto sc      = m_normDist.initialize(randSvc, Rndm::Gauss(0.0, 1.0));
      if (!sc.isSuccess()) {
        return StatusCode::FAILURE;
      }
      // set energy resolution numbers
      for (size_t i = 0; i < u_eRes.size() && i < 3; ++i) {
        eRes[i] = u_eRes[i];
      }

      // using juggler internal units (GeV, mm, radian, ns)
      dyRangeADC = m_dyRangeADC.value() / GeV;
      tRes       = m_tRes.value() / ns;
      stepTDC    = ns / m_resolutionTDC.value();

      // need signal sum
      if (!u_fields.value().empty()) {
        m_geoSvc = service(m_geoSvcName);
        // sanity checks
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

        // get decoders
        try {
          auto id_desc = m_geoSvc->getDetector()->readout(m_readout).idSpec();
          id_mask      = 0;
          std::vector<std::pair<std::string, int>> ref_fields;
          for (size_t i = 0; i < u_fields.size(); ++i) {
            id_mask |= id_desc.field(u_fields[i])->mask();
            // use the provided id number to find ref cell, or use 0
            int ref = i < u_refs.size() ? u_refs[i] : 0;
            ref_fields.emplace_back(u_fields[i], ref);
          }
          ref_mask = id_desc.encode(ref_fields);
          // debug() << fmt::format("Referece id mask for the fields {:#064b}", ref_mask) << endmsg;
        } catch (...) {
          error() << "Failed to load ID decoder for " << m_readout << endmsg;
          return StatusCode::FAILURE;
        }
        id_mask = ~id_mask;
        info() << fmt::format("ID mask in {:s}: {:#064b}", m_readout.value(), id_mask) << endmsg;
        return StatusCode::SUCCESS;
      }

      return StatusCode::SUCCESS;
    }

    StatusCode execute(const EventContext&) const override
    {
      if (!u_fields.value().empty()) {
        signal_sum_digi();
      } else {
        single_hits_digi();
      }
      return StatusCode::SUCCESS;
    }

  private:
    void single_hits_digi() const {
      // input collections
      const auto* const simhits = m_inputHitCollection.get();
      // Create output collections
      auto* rawhits = m_outputHitCollection.createAndPut();
      for (const auto& ahit : *simhits) {
        // Note: juggler internal unit of energy is GeV
        const double eDep    = ahit.getEnergy();

        // apply additional calorimeter noise to corrected energy deposit
        const double eResRel = (eDep > m_threshold)
            ? m_normDist() * std::sqrt(
                  std::pow(eRes[0] / std::sqrt(eDep), 2) +
                  std::pow(eRes[1], 2) +
                  std::pow(eRes[2] / (eDep), 2)
              )
            : 0;

        const double ped    = m_pedMeanADC + m_normDist() * m_pedSigmaADC;
        const long long adc = std::llround(ped +  eDep * (m_corrMeanScale + eResRel) / dyRangeADC * m_capADC);

        double time = std::numeric_limits<double>::max();
        for (const auto& c : ahit.getContributions()) {
          if (c.getTime() <= time) {
            time = c.getTime();
          }
        }
        const long long tdc = std::llround((time + m_normDist() * tRes) * stepTDC);

        rawhits->create(
          ahit.getCellID(),
          (adc > m_capADC.value() ? m_capADC.value() : adc),
          tdc
        );
      }
    }

    void signal_sum_digi() const {
      const auto* const simhits = m_inputHitCollection.get();
      auto* rawhits = m_outputHitCollection.createAndPut();

      // find the hits that belong to the same group (for merging)
      std::unordered_map<long long, std::vector<edm4hep::SimCalorimeterHit>> merge_map;
      for (const auto &ahit : *simhits) {
        int64_t hid = (ahit.getCellID() & id_mask) | ref_mask;
        auto    it  = merge_map.find(hid);

        if (it == merge_map.end()) {
          merge_map[hid] = {ahit};
        } else {
          it->second.push_back(ahit);
        }
      }

      // signal sum
      for (auto &[id, hits] : merge_map) {
        double edep     = hits[0].getEnergy();
        double time     = hits[0].getContributions(0).getTime();
        double max_edep = hits[0].getEnergy();
        // sum energy, take time from the most energetic hit
        for (size_t i = 1; i < hits.size(); ++i) {
          edep += hits[i].getEnergy();
          if (hits[i].getEnergy() > max_edep) {
            max_edep = hits[i].getEnergy();
            for (const auto& c : hits[i].getContributions()) {
              if (c.getTime() <= time) {
                time = c.getTime();
              }
            }
          }
        }

        // safety check
        const double eResRel = (edep > m_threshold)
            ? m_normDist() * eRes[0] / std::sqrt(edep) +
              m_normDist() * eRes[1] +
              m_normDist() * eRes[2] / edep
            : 0;

        double    ped     = m_pedMeanADC + m_normDist() * m_pedSigmaADC;
        unsigned long long adc     = std::llround(ped + edep * (1. + eResRel) / dyRangeADC * m_capADC);
        unsigned long long tdc     = std::llround((time + m_normDist() * tRes) * stepTDC);

        rawhits->create(
          id,
          (adc > m_capADC.value() ? m_capADC.value() : adc),
          tdc
        );
      }
    }
  };
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(CalorimeterHitDigi)

} // namespace Jug::Digi
