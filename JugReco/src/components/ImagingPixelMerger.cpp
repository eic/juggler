/*
 *  A hits converter to prepare dataset for machine learning
 *  It converts hits with (x, y, z, E) to (E, eta, phi) layer by layer
 *  With a defined grid size and ranges, it merge the hits within one grid and drop-off hits
 * out-of-range The capacity of each layer is fixed (padding with zeros), and the hits with least
 * energies that exceed the capacity will be discarded.
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
#include "JugBase/Utilities/UniqueID.hpp"

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

  /** Converter  for ML algorithm input.
   *
   *  A hits converter to prepare dataset for machine learning
   *  It converts hits with (x, y, z, E) to (E, eta, phi) layer by layer
   *  With a defined grid size and ranges, it merge the hits within one grid and drop-off hits
   * out-of-range The capacity of each layer is fixed (padding with zeros), and the hits with least
   * energies that exceed the capacity will be discarded.
   *
   * \ingroup reco
   */
  class ImagingPixelMerger : public GaudiAlgorithm {
  private:
    // Unique identifier for this hit type, based on the algorithm name
    using HitClassificationType = decltype(eic::CalorimeterHitData().type);
    const HitClassificationType m_type;
  public:
    Gaudi::Property<int>                    m_nHits{this, "numberOfHits", 20};
    Gaudi::Property<int>                    m_nLayers{this, "numberOfLayers", 20};
    Gaudi::Property<double>                 m_etaSize{this, "etaSize", 0.001};
    Gaudi::Property<double>                 m_phiSize{this, "phiSize", 0.001};
    DataHandle<eic::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection",
                                                                 Gaudi::DataHandle::Reader, this};
    DataHandle<eic::CalorimeterHitCollection> m_outputHitCollection{"outputHitCollection",
                                                                  Gaudi::DataHandle::Writer, this};

    ImagingPixelMerger(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , m_type{uniqueID<HitClassificationType>(name)}
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
      std::vector<std::unordered_map<std::pair<int, int>, double, pair_hash>> group_hits(m_nLayers);
      for (const auto& h : hits) {
        auto k = h.layer();
        if ((int)k >= m_nLayers) {
          continue;
        }
        double r   = std::sqrt(h.position().x * h.position().x + h.position().y * h.position().y + h.position().z * h.position().z);
        double th  = std::acos(h.position().z / r);
        double eta = -std::log(std::tan(th / 2.));
        double phi = std::atan2(h.position().y, h.position().x);
        // debug() << th << ", " << eta << ", " << phi << endmsg;

        auto& layer = group_hits[k];
        auto  g     = std::pair<int, int>{(eta + 4.) / m_etaSize, (phi + M_PI) / m_phiSize};
        auto  it    = layer.find(g);
        // merge energy
        if (it != layer.end()) {
          it->second += h.energy();
        } else {
          layer[g] = h.energy();
        }
      }

      // convert to data
      struct GridData {
        double energy, eta, phi;
      };
      for (auto [i, layer] : Jug::Utils::Enumerate(group_hits)) {
        std::vector<GridData> grids;
        for (auto& it : layer) {
          double ie = it.first.first;
          double ip = it.first.second;
          grids.push_back(GridData{it.second, ie * m_etaSize - 4., ip * m_phiSize - M_PI});
          // debug() << ie << ", " << ip << ", " << ie*m_etaSize - 4. << ", " << ip*m_phiSize - M_PI
          // << endmsg;
        }
        std::sort(grids.begin(), grids.end(),
                  [](const GridData& g1, const GridData& g2) { return g1.energy < g2.energy; });

        for (size_t k = 0; k < (size_t)m_nHits; ++k) {
          GridData grid{0., 0., 0.};
          if (k < grids.size()) {
            grid = grids[k];
          }
          eic::VectorXYZ pos {eic::VectorPolar(0., 0., grid.phi)};
          // @TODO: This seems incomplete...
          auto h = mhits.create();
          h.ID(k);
          h.layer(i);
          h.type(m_type);
          h.energy(grid.energy);
          h.position(pos);
        }
      }

      return StatusCode::SUCCESS;
    }

  }; // class ImagingPixelMerger

  DECLARE_COMPONENT(ImagingPixelMerger)

} // namespace Jug::Reco

