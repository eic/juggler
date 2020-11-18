#include <algorithm>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/PhysicalConstants.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  class SimpleClustering : public GaudiAlgorithm {
  public:
    using RecHits  = eic::CalorimeterHitCollection;
    using Clusters = eic::ClusterCollection;

    DataHandle<RecHits>     m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<Clusters>    m_outputClusters{"outputClusters", Gaudi::DataHandle::Writer, this};
    Gaudi::Property<double> m_minModuleEdep{this, "minModuleEdep", 5.0*MeV};
    Gaudi::Property<double> m_maxDistance{this, "maxDistance", 20.0*cm};

    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    SimpleClustering(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputHitCollection",   m_inputHitCollection,   "");
        declareProperty("outputClusters",  m_outputClusters,  "Output clusters");
    }

    StatusCode initialize() override
    {
        if (GaudiAlgorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }
        m_geoSvc = service("GeoSvc");
        if (!m_geoSvc) {
            error() << "Unable to locate Geometry Service. "
                    << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
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

      std::vector<eic::CalorimeterHit> hits_A;
      std::vector<eic::CalorimeterHit> hits_B;
      std::vector<eic::CalorimeterHit>& the_hits  = hits_A;
      std::vector<eic::CalorimeterHit>& remaining_hits = hits_B;

      double max_dist = m_maxDistance.value() / mm;

      eic::CalorimeterHit ref_hit;
      ref_hit.energy(0.0);
      for (const auto& h : hits) {
        if (h.energy() > ref_hit.energy()) {
          ref_hit = h;
        }
        the_hits.push_back(h.clone());
      }

      debug() << " max_dist = " << max_dist << endmsg;
      bool continue_clustering = true;

      while (continue_clustering) {

        std::vector<eic::CalorimeterHit> cluster_hits;
        //debug() << " ref_hit" << ref_hit << endmsg;
        // auto ref_hit = *(std::max_element(std::begin(hits), std::end(hits),
        //                                [](const auto& h1, const auto& h2) { return h1.energy() < h2.energy(); }));
        // std::partition_copy(std::begin(hits), std::end(hits), std::begin(first_cluster_hits),
        // std::begin(remaining_hits),
        //                    [=](const auto& h) {
        //                      return (std::hypot(h.x() - ref_hit.x(), h.y() - ref_hit.y(), h.z() - ref_hit.z()) <
        //                              max_dist);
        //                    });

        for (const auto& h : the_hits) {
          if (std::hypot(h.x() - ref_hit.x(), h.y() - ref_hit.y(), h.z() - ref_hit.z()) < max_dist) {
            cluster_hits.push_back(h.clone());
          } else {
            remaining_hits.push_back(h.clone());
          }
        }

        debug() << " cluster size " << cluster_hits.size() << endmsg;
        double total_energy =
            std::accumulate(std::begin(cluster_hits), std::end(cluster_hits), 0.0,
                            [](double t, const eic::CalorimeterHit& h1) { return (t + h1.energy()); });
        debug() << " total_energy = " << total_energy << endmsg;

        eic::Cluster cluster0;
        // debug() << " cluster0 = " << cluster0 << endmsg;
        // eic::Cluster cluster1  = std::accumulate(std::begin(first_cluster_hits), std::end(first_cluster_hits),
        // cluster0,
        for (const auto& h : cluster_hits) {
          cluster0.energy(cluster0.energy() + h.energy());
          cluster0.x(cluster0.x() + h.energy() * h.x() / total_energy);
          cluster0.y(cluster0.y() + h.energy() * h.y() / total_energy);
          cluster0.z(cluster0.z() + h.energy() * h.z() / total_energy);
        }
        //debug() << " cluster = " << cluster0 << endmsg;
        clusters.push_back(cluster0);


        if((remaining_hits.size() > 5) && (clusters.size() < 10) ){
          ref_hit.energy(0.0);
          for (const auto& h : remaining_hits) {
            if (h.energy() > ref_hit.energy()) {
              ref_hit = h;
            }
          }

          std::swap( remaining_hits, the_hits);
          remaining_hits.clear();

        } else {
          continue_clustering = false;
          break;

        }
      }
        debug() << " clusters found: " << clusters.size() << endmsg;

      return StatusCode::SUCCESS;
      }
  };

  DECLARE_COMPONENT(SimpleClustering)

} // namespace Jug::Reco
