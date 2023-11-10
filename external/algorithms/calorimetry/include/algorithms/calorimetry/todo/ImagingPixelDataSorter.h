// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Chao, Sylvester Joosten

/*
 *  A hits-level data sorter to prepare dataset for machine learning
 *
 *  Author: Chao Peng (ANL), 05/04/2021
 */
#include <algorithm>
#include <bitset>
#include <unordered_map>

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
#include "JugBase/Utilities/Utils.hpp"

// Event Model related classes
#include <edm4hep/utils/vector_utils.h>
#include "edm4eic/CalorimeterHit.h"
#include "edm4eic/CalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  /** Hits sorter for ML algorithm input.
   *
   * A hits-level data sorter to prepare dataset for machine learning
   * It sorts the hits by layer and energy with defined sizes (max number of layers and max number of hits per layer)
   * Hits are sorted by energy in a descending order.
   * Out-of-range hits will be discarded and empty slots will be padded with zeros
   *
   * \ingroup reco
   */
  class ImagingPixelDataSorter : public GaudiAlgorithm {
  private:
    Gaudi::Property<int>                        m_nLayers{this, "numberOfLayers", 9};
    Gaudi::Property<int>                        m_nHits{this, "numberOfHits", 50};
    DataHandle<edm4eic::CalorimeterHitCollection>   m_inputHitCollection{"inputHitCollection",
                                                                     Gaudi::DataHandle::Reader, this};
    DataHandle<edm4eic::CalorimeterHitCollection>   m_outputHitCollection{"outputHitCollection",
                                                                      Gaudi::DataHandle::Writer, this};

  public:
    ImagingPixelDataSorter(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& hits = *m_inputHitCollection.get();
      // Create output collections
      auto& mhits = *m_outputHitCollection.createAndPut();

      // group the hits by layer
      std::vector<std::vector<edm4eic::CalorimeterHit>> layer_hits;
      layer_hits.resize(m_nLayers);
      for (const auto& h : hits) {
        auto k = h.getLayer();
        if ((int)k < m_nLayers) {
          layer_hits[k].push_back(h);
        }
      }

      // sort by energy
      for (auto &layer : layer_hits) {
        std::sort(layer.begin(), layer.end(),
          [] (const edm4eic::CalorimeterHit &h1, const edm4eic::CalorimeterHit &h2) {
            return h1.getEnergy() > h2.getEnergy();
          });
      }

      // fill-in the output
      for (size_t k = 0; k < layer_hits.size(); ++k) {
        auto &layer = layer_hits[k];
        for (size_t i = 0; i < (size_t) m_nHits; ++i) {
          // pad zeros if no hits
          if (i >= layer.size()) {
            auto h = mhits.create();
            h.setLayer((int)k);
            h.setEnergy(0.);
          } else {
            auto h = layer[i].clone();
            mhits.push_back(h);
          }
        }
      }

      return StatusCode::SUCCESS;
    }

  }; // class ImagingPixelDataSorter

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(ImagingPixelDataSorter)

} // namespace Jug::Reco
