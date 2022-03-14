
/*
 *  Fix me
 *
 *  Author: Fix me
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
 * It accepts inputs from data sorter that hits are sorted by layers
 * Two different datasets will be combined together following specified rules in handling the layers
 * Supported rules: concatenate, interlayer
 *
 * \ingroup reco
 */
class ImagingPixelDataShaper: public GaudiAlgorithm {
public:
  DataHandle<eicd::CalorimeterHitCollection> m_inputHits{"inputHits", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::CalorimeterHitCollection> m_outputHits{"outputHits", Gaudi::DataHandle::Writer, this};

  ImagingPixelDataShaper(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHits", m_inputHits, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto hits = m_inputHits.get();

    float output[20][29][5];  // fill array with hits

    return StatusCode::SUCCESS;
  }

}; // class ImagingPixelDataShaper

DECLARE_COMPONENT(ImagingPixelDataShaper)

} // namespace Jug::Reco

