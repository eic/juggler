/*
 *  A hits merger for ecal barrel to prepare dataset for machine learning
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
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "eicd/VectorPolar.h"
#include "eicd/VectorXYZ.h"
#include "eicd/CalorimeterHit.h"
#include "eicd/CalorimeterHitCollection.h"

using namespace Gaudi::Units;

struct pair_hash {
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
  class ImagingPixelMerger : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:
    Gaudi::Property<float>                  m_etaSize{this, "etaSize", 0.001};
    Gaudi::Property<float>                  m_phiSize{this, "phiSize", 0.001};
    DataHandle<eic::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection",
                                                                 Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection> m_outputHitCollection{"outputHitCollection",
                                                                  Gaudi::DataHandle::Writer, this};

    ImagingPixelMerger(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin<>(name, info())
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

      // @TODO: add timing information
      // group the hits by grid per layer
      struct GridData {
        float rc, energy;
      };
      // @TODO: remove this hard-coded value
      int max_nlayers = 50;
      std::vector<std::unordered_map<std::pair<int, int>, GridData, pair_hash>> group_hits(max_nlayers);
      for (const auto& h : hits) {
        auto k = h.layer();
        if ((int)k >= max_nlayers) {
          continue;
        }
        auto& layer = group_hits[k];
        auto& pos = h.position();

        // cylindrical r
        float rc  = std::sqrt(pos.x * pos.x + pos.y * pos.y);
        float eta = pos.eta();
        float phi = pos.phi();

        auto  grid  = std::pair<int, int>{pos2grid(eta, m_etaSize), pos2grid(phi, m_phiSize)};
        auto  it    = layer.find(grid);
        // merge energy
        if (it != layer.end()) {
          it->second.energy += h.energy();
        } else {
          layer[grid] = GridData{rc, h.energy()};
        }
      }

      // convert grid data back to hits
      int k = 0;
      for (auto [i, layer] : Jug::Utils::Enumerate(group_hits)) {
        for (auto [grid, data] : layer) {
          float eta = grid2pos(grid.first, m_etaSize);
          float phi = grid2pos(grid.second, m_phiSize);
          float theta = std::atan(std::exp(-eta))*2.;
          double z = cotan(theta)*data.rc;
          float r = std::sqrt(data.rc*data.rc + z*z);

          eic::VectorXYZ pos {eic::VectorPolar(r, theta, phi)};
          auto h = mhits.create();
          h.ID({k++, algorithmID()});
          h.layer(i);
          h.energy(data.energy);
          h.position(pos);
        }
      }
      return StatusCode::SUCCESS;
    }

  private:
    inline int pos2grid(float pos, float size, float offset = 0.) { return std::floor((pos - offset)/size); }
    inline int grid2pos(float grid, float size, float offset = 0.) { return (grid + 0.5)*size + offset; }
    inline float cotan(float theta) { if (std::abs(std::sin(theta)) < 1e-6) { return 0.; } else { return 1./std::tan(theta); } }
  }; // class ImagingPixelMerger

  DECLARE_COMPONENT(ImagingPixelMerger)

} // namespace Jug::Reco

