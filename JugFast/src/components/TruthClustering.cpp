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
class TruthClustering : public GaudiAlgorithm, AlgorithmIDMixin<> {
public:
  using RecHits       = eic::CalorimeterHitCollection;
  using ProtoClusters = eic::ProtoClusterCollection;
  using TruthHits     = edm4hep::SimCalorimeterHitCollection;

  DataHandle<RecHits> m_inputHits{"inputHits", Gaudi::DataHandle::Reader, this};
  DataHandle<TruthHits> m_mcHits{"mcHits", Gaudi::DataHandle::Reader, this};
  DataHandle<ProtoClusters> m_outputProtoClusters{"outputProtoClusters", Gaudi::DataHandle::Writer, this};

  // Monte Carlo particle source identifier
  const int32_t kMonteCarloSource{uniqueID<int32_t>("MCParticles")};

  TruthClustering(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin<>(name, info()) {
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

    std::vector<std::pair<uint32_t, eic::ConstCalorimeterHit>> the_hits;
    std::vector<std::pair<uint32_t, eic::ConstCalorimeterHit>> remaining_hits;

    // Map mc track ID to protoCluster index
    std::map<int32_t, int32_t> protoIndex;

    // Loop over al calorimeter hits and sort per mcparticle
    int32_t clusterID = 0;
    for (uint32_t i = 0; i < hits.size(); ++i) {
      const auto& hit       = hits[i];
      const auto& mcHit     = mc[hit.ID().value];
      const int32_t trackID = mcHit.truth().trackID;
      // Create a new protocluster if we don't have one for this trackID
      if (!protoIndex.count(trackID)) {
        auto pcl = proto.create();
        pcl.ID({clusterID++, algorithmID()});
        protoIndex[trackID] = proto.size() - 1;
      }
      // Add hit to the appropriate protocluster
      proto[protoIndex[trackID]].addhits({hit.ID(), i, 1.});
    }
    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT(TruthClustering)

} // namespace Jug::Fast
