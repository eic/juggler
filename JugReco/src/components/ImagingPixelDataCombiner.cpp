/*
 *  A hits-level data combiner to combine two datasets into one for machine learning
 *
 *  Author: Chao Peng (ANL), 05/04/2021
 */
#include <algorithm>
#include <bitset>
#include <fmt/format.h>
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
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/vector_utils.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

/** Hits combiner for ML algorithm input.
 *
 * A hits-level data combiner to combine two datasets into one for machine learning
 * It accepts inputs from data x that hits are sorted by layers
 * Two different datasets will be combined together following specified rules in handling the layers
 * Supported rules: concatenate, interlayer
 *
 * \ingroup reco
 */
class ImagingPixelDataCombiner : public GaudiAlgorithm {
private:
  Gaudi::Property<int> m_layerIncrement{this, "layerIncrement", 0};
  Gaudi::Property<std::string> m_rule{this, "rule", "concatenate"};
  DataHandle<eicd::CalorimeterHitCollection> m_inputHits1{"inputHits1", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::CalorimeterHitCollection> m_inputHits2{"inputHits2", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::CalorimeterHitCollection> m_outputHits{"outputHits", Gaudi::DataHandle::Writer, this};
  std::vector<std::string> supported_rules{"concatenate", "interlayer"};

public:
  ImagingPixelDataCombiner(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHits1", m_inputHits1, "");
    declareProperty("inputHits2", m_inputHits2, "");
    declareProperty("outputHits", m_outputHits, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    if (std::find(supported_rules.begin(), supported_rules.end(), m_rule.value()) == supported_rules.end()) {
      error() << fmt::format("unsupported rule: {}, please choose one from [{}]", m_rule.value(),
                             fmt::join(supported_rules, ", "))
              << endmsg;
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto* const hits1 = m_inputHits1.get();
    const auto* const hits2 = m_inputHits2.get();
    std::vector<const eicd::CalorimeterHitCollection*> inputs{hits1, hits2};
    // Create output collections
    auto* mhits = m_outputHits.createAndPut();

    // concatenate
    if (m_rule.value() == supported_rules[0]) {
      for (int i = 0; i < (int)inputs.size(); ++i) {                                   
        auto coll = inputs[i];
        for (auto hit : *coll) {
          eicd::CalorimeterHit h2{
              hit.getCellID(),    hit.getEnergy(),   hit.getEnergyError(), hit.getTime(),
              hit.getTimeError(), hit.getPosition(), hit.getDimension(),   hit.getLayer() + m_layerIncrement * i,
              hit.getSector(),    hit.getLocal(),
          };
          mhits->push_back(h2);
        }
      }
      // interlayer
      // @NOTE: it assumes the input hits are sorted by layers
    } else if (m_rule.value() == supported_rules[1]) {
      std::vector<int> indices{0, 0};
      int curr_coll   = 0;
      bool init_layer = false;
      int curr_layer  = 0;
      // int curr_ihit = 0;
      while (indices[0] < (int)hits1->size() || indices[1] < (int)hits2->size()) {
        // cyclic index
        if (curr_coll >= (int)inputs.size()) {
          curr_coll -= (int)inputs.size();
        }

        // merge hits
        int& i                 = indices[curr_coll];
        const auto* const coll = inputs[curr_coll];

        // reach this collection's end
        if (i >= (int)coll->size()) {
          curr_coll++;
          init_layer = false;
          // curr_ihit = 0;
          // info() << "collection end" << endmsg;
          continue;
        }

        auto hit = (*coll)[i];
        if (!init_layer) {
          curr_layer = hit.getLayer();
          init_layer = true;
        }

        // reach this layer's end
        if (curr_layer != hit.getLayer()) {
          curr_coll++;
          init_layer = false;
          // curr_ihit = 0;
          // info() << "layer end : " << curr_layer << " != " << hit.getLayer() << endmsg;
          continue;
        }

        // push hit, increment of index
        eicd::CalorimeterHit h2{
            hit.getCellID(),    hit.getEnergy(),   hit.getEnergyError(), hit.getTime(),
            hit.getTimeError(), hit.getPosition(), hit.getDimension(),   hit.getLayer() + m_layerIncrement * curr_coll,
            hit.getSector(),    hit.getLocal()};
        mhits->push_back(h2);
        i++;
        // info() << curr_coll << ": " << curr_ihit ++ << endmsg;
      }
      // info() << hits1->size() << ", " << hits2->size() << endmsg;
    }

    return StatusCode::SUCCESS;
  }

}; // class ImagingPixelDataCombiner

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ImagingPixelDataCombiner)

} // namespace Jug::Reco

