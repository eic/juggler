
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
    const auto& hits = *m_inputHits.get();  // 1d

    // group the hits by layer
    std::vector<std::vector<eicd::CalorimeterHit>> layer_hits;  // create a vector of vector of hits
    layer_hits.resize(int m_nLayers = 29);  // create 29 layers/rows 
    for (const auto& h : hits) {
      auto k = h.getLayer();
      if ((int)k < m_nLayers) {
        layer_hits[k].push_back(h);  // take each hit and append it to vector above 
      }
    }

    // sort by energy. Descending energies for each layer. 
    for (auto &layer : layer_hits) {
      std::sort(layer.begin(), layer.end(),
        [] (const eicd::CalorimeterHit &h1, const eicd::CalorimeterHit &h2) {
          return h1.getEnergy() > h2.getEnergy();
        });
      layer.resize(int topHits = 20);  // get only the 20 most energetic hits for each layer
    }

    int layers = 29;
    int nHits = 20;
    int features = 5;
    
    float output[layer_hits[0].size()][layer_hits[1].size] = {};

    for (int i = 0; i < layer_hits.size(); i++) {  // 29 layers
      for (int j = 0; j < layer_hits[0].size(); j++) {  // 20 most energetic hits
          x = layer_hits[i][j].getPosition().x;
          y = layer_hits[i][j].getPosition().y;
          z = layer_hits[i][j].getPosition().z;

          rc = sqrt(pow(x, 2) + pow(y, 2)); 
          r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
          theta = acos(z/r);
          phi = atan2(y, x)
          eta = -1*log10(tan(theta/2))

          layer_type = floor(layer_hits[i][j].getLayer());
          energy = layer_hits[i][j].getEnergy();
             
          output[i][j] = {layer_type, energy, rc, eta, phi}
        }
      }
    }
      

    return StatusCode::SUCCESS;
  }

}; // class ImagingPixelDataShaper

DECLARE_COMPONENT(ImagingPixelDataShaper)

} // namespace Jug::Reco

