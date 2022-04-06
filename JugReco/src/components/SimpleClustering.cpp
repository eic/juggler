// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck, Chao
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

// Event Model related classes
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ProtoClusterCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/vector_utils.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  /** Simple clustering algorithm.
   *
   * \ingroup reco
   */
  class SimpleClustering : public GaudiAlgorithm {
  private:
    using RecHits  = eicd::CalorimeterHitCollection;
    using ProtoClusters = eicd::ProtoClusterCollection;
    using Clusters = eicd::ClusterCollection;

    DataHandle<RecHits>       m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<ProtoClusters> m_outputProtoClusters{"outputProtoCluster", Gaudi::DataHandle::Writer, this};
    DataHandle<Clusters>      m_outputClusters{"outputClusterCollection", Gaudi::DataHandle::Writer, this};

    Gaudi::Property<std::string> m_mcHits{this, "mcHits", ""};

    Gaudi::Property<double>   m_minModuleEdep{this, "minModuleEdep", 5.0 * MeV};
    Gaudi::Property<double>   m_maxDistance{this, "maxDistance", 20.0 * cm};

    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    // Optional handle to MC hits
    std::unique_ptr<DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputMC;

  public:
    SimpleClustering(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc) {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputProtoClusterCollection", m_outputClusters, "Output proto clusters");
      declareProperty("outputClusterCollection", m_outputClusters, "Output clusters");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      // Initialize the MC input hit collection if requested
      if (m_mcHits != "") {
        m_inputMC =
            std::make_unique<DataHandle<edm4hep::SimCalorimeterHitCollection>>(m_mcHits, Gaudi::DataHandle::Reader, this);
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

      // Optional MC data
      // FIXME no connection between cluster and truth in edm4hep
      //const edm4hep::SimCalorimeterHitCollection* mcHits = nullptr;
      //if (m_inputMC) {
      //  mcHits = m_inputMC->get();
      //}

      std::vector<std::pair<uint32_t, eicd::CalorimeterHit>> the_hits;
      std::vector<std::pair<uint32_t, eicd::CalorimeterHit>> remaining_hits;

      double max_dist   = m_maxDistance.value() / mm;
      double min_energy = m_minModuleEdep.value() / GeV;

      eicd::CalorimeterHit ref_hit;
      bool have_ref = false;
      // Collect all our hits, and get the highest energy hit
      {
        uint32_t idx  = 0;
        for (const auto& h : hits) {
          if (!have_ref || h.getEnergy() > ref_hit.getEnergy()) {
            ref_hit  = h;
            have_ref = true;
          }
          the_hits.emplace_back(idx, h);
          idx += 1;
        }
      }

      if (msgLevel(MSG::DEBUG)) {
        debug() << " max_dist = " << max_dist << endmsg;
      }

      while (have_ref && ref_hit.getEnergy() > min_energy) {

        std::vector<std::pair<uint32_t, eicd::CalorimeterHit>> cluster_hits;

        for (const auto& [idx, h] : the_hits) {
          if (eicd::magnitude(h.getPosition() - ref_hit.getPosition()) < max_dist) {
            cluster_hits.emplace_back(idx, h);
          } else {
            remaining_hits.emplace_back(idx, h);
          }
        }

        double total_energy = std::accumulate(
            std::begin(cluster_hits), std::end(cluster_hits), 0.0,
            [](double t, const std::pair<uint32_t, eicd::CalorimeterHit>& h1) { return (t + h1.second.getEnergy()); });

        if (msgLevel(MSG::DEBUG)) {
          debug() << " total_energy = " << total_energy << endmsg;
          debug() << " cluster size " << cluster_hits.size() << endmsg;
        }
        auto cl = clusters.create();
        cl.setNhits(cluster_hits.size());
        auto pcl = proto.create();
        for (const auto& [idx, h] : cluster_hits) {
          cl.setEnergy(cl.getEnergy() + h.getEnergy());
          cl.setPosition(cl.getPosition() + (h.getPosition() * h.getEnergy() / total_energy));
          pcl.addToHits(h);
          pcl.addToWeights(1);
        }
        // Optionally store the MC truth associated with the first hit in this cluster
        // FIXME no connection between cluster and truth in edm4hep
        //if (mcHits) {
        //  const auto& mc_hit = (*mcHits)[ref_hit.ID().value];
        //  cl.mcID({mc_hit.truth().trackID, m_kMonteCarloSource});
        //}

        have_ref = false;
        if ((remaining_hits.size() > 5) && (clusters.size() < 10)) {
          for (const auto& [idx, h] : remaining_hits) {
            if (!have_ref || h.getEnergy() > ref_hit.getEnergy()) {
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

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  DECLARE_COMPONENT(SimpleClustering)

} // namespace Jug::Reco
