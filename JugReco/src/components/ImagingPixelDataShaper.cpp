
/*
 *  Fix me
 *
 *  Author: Fix me
 */
#include <algorithm>
#include <bitset>
#include <fmt/format.h>
#include <unordered_map>
#include <cmath>
#include <iostream>

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
using namespace std;

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
    
    auto coll = hits;
    int count = 1;
    for (auto hit : *coll) {
      cout << count << endl;
      // cout << hit.getLayer().size_t << endl;
      count++;
      // eicd::CalorimeterHit h2{
      //   hit.getEnergy(),  hit.getPosition(), hit.getLayer(),
      //   };
    }

    

    // int layers = 29;
    // int nHits = 20;
    // int features = 5;
    // float output[layers][nHits][features];

    // for (int i = 0; i < layers; i++) {  // 29 layers
    //   for (int j = 0; j < nHits; j++) {  // 20 most energetic hits
    //     for (int k = 0; k < features; k++){  // 5 features of each hit
    //       x = hits.getPosition().x;
    //       y = hits.getPosition().y;
    //       z = hits.getPosition().z;

    //       rc = sqrt(pow(x, 2) + pow(y, 2));
    //       r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    //       theta = acos(z/r);
    //       phi = atan2(hits.getPosition.y, hits.getPosition.x)
    //       eta = -1*log10(tan(theta/ 2.))

    //       layer_type = floor(hits.getLayer());
    //       energy = hits.getEnergy();
    //       output[i][j][k] = {layer_type, energy, rc, eta, phi}
    //     }
    //   }
    // }
      

    return StatusCode::SUCCESS;
  }

}; // class ImagingPixelDataShaper

DECLARE_COMPONENT(ImagingPixelDataShaper)

} // namespace Jug::Reco

