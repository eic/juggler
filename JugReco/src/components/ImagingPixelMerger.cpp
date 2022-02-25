/*
 *  A hits merger for ecal barrel to prepare dataset for machine learning
 *
 *  Author: Chao Peng (ANL), 05/04/2021
 *
 *  SJJ: @TODO: this should really be a clustering algorithm, as it doesn't give us
 *              fully consistent hits anymore (e.g. cellID, local position) (meaning it's 
 *              not generically useful as a reconstruction algorithm.
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
#include "eicd/CalorimeterHitCollection.h"
#include <eicd/vector_utils.h>

using namespace Gaudi::Units;

struct PairHashFunction {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2>& pair) const
  {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};


namespace Jug::Reco {

  /** Hits merger for ML algorithm input.
   *
   * A hits merger to prepare dataset for machine learning
   * It merges hits with a defined grid size.
   * Merged hits will be relocated to the grid center and the energies will be summed.
   *
   * \ingroup reco
   */
  class ImagingPixelMerger : public GaudiAlgorithm {
  public:
    Gaudi::Property<float>                  m_etaSize{this, "etaSize", 0.001};
    Gaudi::Property<float>                  m_phiSize{this, "phiSize", 0.001};
    DataHandle<eicd::CalorimeterHitCollection> m_inputHits{"inputHits", Gaudi::DataHandle::Reader, this};
    DataHandle<eicd::CalorimeterHitCollection> m_outputHits{"outputHits", Gaudi::DataHandle::Writer, this};

    ImagingPixelMerger(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHits", m_inputHits, "");
      declareProperty("outputHits", m_outputHits, "");
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
      const auto& hits = *m_inputHits.get();
      // Create output collections
      auto& ohits = *m_outputHits.createAndPut();

      // @TODO: add timing information
      // group the hits by grid per layer
      struct GridData {
        unsigned nHits;
        float rc;
        float energy;
        float energyError;
        float time;
        float timeError;
        int sector; // sector associated with one of the merged hits
      };
      // @TODO: remove this hard-coded value
      int max_nlayers = 50;
      std::vector<std::unordered_map<std::pair<int, int>, GridData, PairHashFunction>> group_hits(max_nlayers);
      for (const auto& h : hits) {
        auto k = h.layer();
        if ((int)k >= max_nlayers) {
          continue;
        }
        auto& layer = group_hits[k];
        const auto& pos = h.position();

        // cylindrical r
        const float rc  = eicd::magnitudeTransverse(pos);
        const double eta = eicd::eta(pos);
        const double phi = eicd::angleAzimuthal(pos);

        const auto grid = std::pair<int, int>{pos2grid(eta, m_etaSize), pos2grid(phi, m_phiSize)};
        auto it         = layer.find(grid);
        // merge energy
        if (it != layer.end()) {
          auto& data = it->second;
          data.nHits += 1;
          data.energy += h.energy();
          data.energyError += h.energyError() * h.energyError();
          data.time += h.time();
          data.timeError += h.timeError() * h.timeError();
        } else {
          layer[grid] = GridData{
              1,         rc, h.energy(), h.energyError() * h.energyError(), h.time(), h.timeError() * h.timeError(),
              h.sector()};
        }
      }

      // convert grid data back to hits
      for (const auto& [i, layer] : Jug::Utils::Enumerate(group_hits)) {
        for (const auto& [grid, data] : layer) {
          const double eta = grid2pos(grid.first, m_etaSize);
          const double phi = grid2pos(grid.second, m_phiSize);
          const double theta = eicd::etaToAngle(eta);
          const double z   = cotan(theta) * data.rc;
          const float r    = std::hypot(data.rc, z);
          const auto pos = eicd::sphericalToVector(r, theta, phi);
          auto oh = ohits.create();
          oh.energy(data.energy);
          oh.energyError(std::sqrt(data.energyError));
          oh.time(data.time / data.nHits);
          oh.timeError(std::sqrt(data.timeError));
          oh.position(pos);
          oh.layer(i);
          oh.sector(data.sector);
        }
      }
      return StatusCode::SUCCESS;
    }

  private:
    int pos2grid(float pos, float size, float offset = 0.) const { return std::floor((pos - offset) / size); }
    int grid2pos(float grid, float size, float offset = 0.) const { return (grid + 0.5) * size + offset; }
    float cotan(float theta) const {
      if (std::abs(std::sin(theta)) < 1e-6) {
        return 0.;
      } else {
        return 1. / std::tan(theta);
      }
    }
  }; // class ImagingPixelMerger

  DECLARE_COMPONENT(ImagingPixelMerger)

} // namespace Jug::Reco

