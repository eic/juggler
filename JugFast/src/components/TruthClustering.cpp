// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Whitney Armstrong, Wouter Deconinck

#include <algorithm>

#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>

// Event Model related classes
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/ProtoClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Fast {

/** Truth clustering algorithm.
 *
 * \ingroup reco
 */
class TruthClustering : public Gaudi::Algorithm {
private:
  mutable DataHandle<edm4eic::CalorimeterHitCollection> m_inputHits{"inputHits", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::SimCalorimeterHitCollection> m_mcHits{"mcHits", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4eic::ProtoClusterCollection> m_outputProtoClusters{"outputProtoClusters", Gaudi::DataHandle::Writer, this};

public:
  TruthClustering(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {
    declareProperty("inputHits", m_inputHits, "Input calorimeter reco hits");
    declareProperty("mcHits", m_mcHits, "Input truth hits");
    declareProperty("outputProtoClusters", m_outputProtoClusters, "Output proto clusters");
  }

  StatusCode initialize() override {
    if (Gaudi::Algorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    // input collections
    const auto& hits = *m_inputHits.get();
    const auto& mc   = *m_mcHits.get();
    // Create output collections
    auto& proto = *m_outputProtoClusters.createAndPut();

    // Map mc track ID to protoCluster index
    std::map<int32_t, int32_t> protoIndex;

    // Loop over al calorimeter hits and sort per mcparticle
    for (const auto& hit : hits) {
      const auto& mcHit     = mc[hit.getObjectID().index];
      const auto& trackID   = mcHit.getContributions(0).getParticle().getObjectID().index;
      // Create a new protocluster if we don't have one for this trackID
      if (protoIndex.count(trackID) == 0) {
        auto pcl = proto.create();
        protoIndex[trackID] = proto.size() - 1;
      }
      // Add hit to the appropriate protocluster
      proto[protoIndex[trackID]].addToHits(hit);
      proto[protoIndex[trackID]].addToWeights(1);
    }
    return StatusCode::SUCCESS;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(TruthClustering)

} // namespace Jug::Fast
