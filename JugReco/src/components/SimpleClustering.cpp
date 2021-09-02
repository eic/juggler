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

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ProtoClusterCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  /** Simple clustering algorithm.
   *
   * \ingroup reco
   */
  class SimpleClustering : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:
    using RecHits  = eic::CalorimeterHitCollection;
    using ProtoClusters = eic::ProtoClusterCollection;
    using Clusters = eic::ClusterCollection;

    DataHandle<RecHits>       m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<ProtoClusters> m_outputProtoClusters{"outputProtoCluster", Gaudi::DataHandle::Writer, this};
    DataHandle<Clusters>      m_outputClusters{"outputClusterCollection", Gaudi::DataHandle::Writer, this};

    Gaudi::Property<double>   m_minModuleEdep{this, "minModuleEdep", 5.0 * MeV};
    Gaudi::Property<double>   m_maxDistance{this, "maxDistance", 20.0 * cm};

    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    SimpleClustering(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin<>(name, info()) {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputProtoClusterCollection", m_outputClusters, "Output proto clusters");
      declareProperty("outputClusterCollection", m_outputClusters, "Output clusters");
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
      auto& proto = *m_outputProtoClusters.createAndPut();
      auto& clusters = *m_outputClusters.createAndPut();

      std::vector<eic::ConstCalorimeterHit> the_hits;
      std::vector<eic::ConstCalorimeterHit> remaining_hits;

      double max_dist   = m_maxDistance.value() / mm;
      double min_energy = m_minModuleEdep.value() / GeV;

      eic::ConstCalorimeterHit ref_hit;
      bool                have_ref = false;
      for (const auto& h : hits) {
        //const eic::CalorimeterHit h = ch.clone();
        if (!have_ref || h.energy() > ref_hit.energy()) {
          ref_hit  = h;
          have_ref = true;
        }
        the_hits.push_back(h);
      }

      if (msgLevel(MSG::DEBUG)) {
        debug() << " max_dist = " << max_dist << endmsg;
      }

      while (have_ref && ref_hit.energy() > min_energy) {

        std::vector<eic::ConstCalorimeterHit> cluster_hits;

        for (const auto& h : the_hits) {
          if (std::hypot(h.position().x - ref_hit.position().x, h.position().y - ref_hit.position().y,
                         h.position().z - ref_hit.position().z) < max_dist) {
            cluster_hits.push_back(h);
          } else {
            remaining_hits.push_back(h); }
        }

        double total_energy =
            std::accumulate(std::begin(cluster_hits), std::end(cluster_hits), 0.0,
                            [](double t, const eic::ConstCalorimeterHit& h1) { return (t + h1.energy()); });

        if (msgLevel(MSG::DEBUG)) {
          debug() << " total_energy = " << total_energy << endmsg;
          debug() << " cluster size " << cluster_hits.size() << endmsg;
        }

        eic::Cluster cl;
        cl.ID({clusters.size(), algorithmID()});
        cl.nhits(cluster_hits.size());
        for (const auto& h : cluster_hits) {
          cl.energy(cl.energy() + h.energy());
          cl.position(cl.position().add(h.position().scale(h.energy()/total_energy)));
          proto.create(h.ID(), cl.ID());
        }
        clusters.push_back(cl);

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
      if (msgLevel(MSG::DEBUG)) {
        debug() << " clusters found: " << clusters.size() << endmsg;
      }

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(SimpleClustering)

} // namespace Jug::Reco
