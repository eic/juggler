#include <algorithm>

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

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  class SimpleClustering : public GaudiAlgorithm {
  public:
    using RecHits  = eic::CalorimeterHitCollection;
    using Clusters = eic::ClusterCollection;

    DataHandle<RecHits> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<Clusters>    m_outputClusters{"outputClusters", Gaudi::DataHandle::Writer, this};
    Gaudi::Property<double> m_minModuleEdep{this, "minModuleEdep", 5.0 * MeV};
    Gaudi::Property<double> m_maxDistance{this, "maxDistance", 20.0 * cm};

    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    SimpleClustering(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputClusters", m_outputClusters, "Output clusters");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& hits = *m_inputHitCollection.get();
      // Create output collections
      auto& clusters = *m_outputClusters.createAndPut();

      std::vector<eic::CalorimeterHit> the_hits;
      std::vector<eic::CalorimeterHit> remaining_hits;

      double max_dist   = m_maxDistance.value() / mm;
      double min_energy = m_minModuleEdep.value() / GeV;

      eic::CalorimeterHit ref_hit;
      bool                have_ref = false;
      for (const auto& ch : hits) {
        const eic::CalorimeterHit h{
            ch.cellID(), ch.clusterID(), ch.layerID(), ch.sectorID(),  ch.energy(),
            ch.time(),   ch.position(),  ch.local(),   ch.dimension(), 1};
        if (!have_ref || h.energy() > ref_hit.energy()) {
          ref_hit  = h;
          have_ref = true;
        }
        the_hits.push_back(h);
      }

      debug() << " max_dist = " << max_dist << endmsg;

      while (have_ref && ref_hit.energy() > min_energy) {

        std::vector<eic::CalorimeterHit> cluster_hits;

        for (const auto& h : the_hits) {
          if (std::hypot(h.x() - ref_hit.x(), h.y() - ref_hit.y(), h.z() - ref_hit.z()) <
              max_dist) {
            cluster_hits.push_back(h);
          } else {
            remaining_hits.push_back(h);
          }
        }

        debug() << " cluster size " << cluster_hits.size() << endmsg;
        double total_energy = std::accumulate(
            std::begin(cluster_hits), std::end(cluster_hits), 0.0,
            [](double t, const eic::CalorimeterHit& h1) { return (t + h1.energy()); });
        debug() << " total_energy = " << total_energy << endmsg;

        eic::Cluster cluster0;
        for (const auto& h : cluster_hits) {
          cluster0.energy(cluster0.energy() + h.energy());
          cluster0.x(cluster0.x() + h.energy() * h.x() / total_energy);
          cluster0.y(cluster0.y() + h.energy() * h.y() / total_energy);
          cluster0.z(cluster0.z() + h.energy() * h.z() / total_energy);
        }
        clusters.push_back(cluster0);

        have_ref = false;
        if ((remaining_hits.size() > 5) && (clusters.size() < 10)) {
          for (const auto& h : remaining_hits) {
            if (!have_ref || h.energy() > ref_hit.energy()) {
              ref_hit  = h;
              have_ref = true;
            }
          }

          std::swap(remaining_hits, the_hits);
          remaining_hits.clear();
        }
      }
      debug() << " clusters found: " << clusters.size() << endmsg;

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(SimpleClustering)

} // namespace Jug::Reco
