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
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ProtoClusterCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Fast {

/** Truth clustering algorithm.
 *
 * \ingroup reco
 */
class TruthClustering : public GaudiAlgorithm {
private:
  DataHandle<eicd::CalorimeterHitCollection> m_inputHits{"inputHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection> m_mcHits{"mcHits", Gaudi::DataHandle::Reader, this};
  DataHandle<eicd::ProtoClusterCollection> m_outputProtoClusters{"outputProtoClusters", Gaudi::DataHandle::Writer, this};

public:
  TruthClustering(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHits", m_inputHits, "Input calorimeter reco hits");
    declareProperty("mcHits", m_mcHits, "Input truth hits");
    declareProperty("outputProtoClusters", m_outputProtoClusters, "Output proto clusters");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
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
      const auto& trackID   = mcHit.getContributions(0).getParticle().id();
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
