#include <limits>
#include <numbers>

#include <fmt/format.h>
// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "eicd/ClusterCollection.h"
#include "eicd/MergedClusterRelationsCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

/** Simple algorithm to merge the energy measurement from cluster1 with the position
 * measurement of cluster2 (in case matching clusters are found). If not, it will
 * propagate the raw cluster from cluster1 or cluster2
 *
 * Matching occurs based on the cluster phi, z and E variables, with tolerances
 * defined in the options file. A negative tolerance effectively disables
 * a check. The energy tolerance is defined as a relative number (e.g. .1)
 *
 * In case of ambiguity the closest cluster is merged.
 *
 * \ingroup reco
 */
class EnergyPositionClusterMerger : public GaudiAlgorithm, public AlgorithmIDMixin<> {
public:
  // Input
  DataHandle<eic::ClusterCollection> m_energyClusters{"energyClusters", Gaudi::DataHandle::Reader, this};
  DataHandle<eic::ClusterCollection> m_positionClusters{"positionClusters", Gaudi::DataHandle::Reader, this};
  // Output
  DataHandle<eic::ClusterCollection> m_outputClusters{"outputClusters", Gaudi::DataHandle::Writer, this};
  DataHandle<eic::MergedClusterRelationsCollection> m_relations{"outputCluster", Gaudi::DataHandle::Writer,
                                                                this}; // namespace Jug::Reco
  // Negative values mean the tolerance check is disabled
  Gaudi::Property<double> m_zToleranceUnits{this, "zTolerance", -1 * cm};
  Gaudi::Property<double> m_phiToleranceUnits{this, "phiTolerance", 20 * degree};
  Gaudi::Property<double> m_energyRelTolerance{this, "energyRelTolerance", 0.3};
  // Unitless (GeV/mm/ns/rad) versions of these tolerances
  double m_zTolerance;
  double m_phiTolerance;

public:
  EnergyPositionClusterMerger(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info()) {
    declareProperty("energyClusters", m_energyClusters, "Cluster collection with good energy precision");
    declareProperty("positionClusters", m_positionClusters, "Cluster collection with good position precision");
    declareProperty("outputClusters", m_outputClusters, "");
    declareProperty("outputCluster", m_relations, "");
  }

  StatusCode initialize() override {
    m_zTolerance   = m_zToleranceUnits / mm;
    m_phiTolerance = m_phiTolerance / rad;
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input
    const auto& e_clus   = *(m_energyClusters.get());
    const auto& pos_clus = *(m_positionClusters.get());
    // output
    auto& merged    = *(m_outputClusters.createAndPut());
    auto& relations = *(m_relations.createAndPut());

    std::vector<bool> consumed(e_clus.size(), false);

    int idx = 0; // merged cluster ID
    // use position clusters as starting point
    for (const auto& pc : pos_clus) {
      // check if we find a good match
      int best_match    = -1;
      double best_delta = std::numeric_limits<double>::max();
      for (size_t ie = 0; ie < e_clus.size(); ++ie) {
        if (consumed[ie]) {
          continue;
        }
        const auto& ec = e_clus[ie];
        // 1. stop if not within tolerance
        //    (make sure to handle rollover of phi properly)
        double dphi = pc.position().phi() - ec.position().phi();
        if (std::abs(dphi) > M_PI) {
          dphi = std::abs(dphi) - M_PI;
        }
        if ((m_energyRelTolerance > 0 &&
             (ec.energy() <= 0 || fabs((pc.energy() - ec.energy()) / ec.energy()) > m_energyRelTolerance)) ||
            (m_zTolerance > 0 && fabs(pc.position().z - ec.position().z) > m_zTolerance) ||
            (m_phiTolerance > 0 && dphi > m_phiTolerance)) {
          continue;
        }
        // --> if we get here, we have a match within tolerance. Now treat the case
        //     where we have multiple matches. In this case take the one with the closest
        //     energies.
        // 2. best match?
        const double delta = fabs(pc.energy() - ec.energy());
        if (delta < best_delta) {
          best_delta = delta;
          best_match = ie;
        }
      }
      // Create a merged cluster if we find a good match
      if (best_match >= 0) {
        const auto& ec = e_clus[best_match];
        auto new_clus  = merged.create();
        new_clus.ID({idx++, algorithmID()});
        new_clus.energy(ec.energy());
        new_clus.energyError(ec.energyError());
        new_clus.time(pc.time());
        new_clus.nhits(pc.nhits() + ec.nhits());
        new_clus.position(pc.position());
        new_clus.positionError(pc.positionError());
        new_clus.radius(pc.radius());
        new_clus.skewness(pc.skewness());
        auto rel = relations.create();
        rel.clusterID(new_clus.ID());
        rel.size(2);
        rel.parent()[0] = pc.ID();
        rel.parent()[1] = ec.ID();
        // label our energy cluster as consumed
        consumed[best_match] = true;
        if (msgLevel(MSG::DEBUG)) {
          debug() << fmt::format("Matched position cluster {} with energy cluster {}\n", pc.ID().value, ec.ID().value) << endmsg;
          debug() << fmt::format("  - Position cluster: (E: {}, phi: {}, z: {})", pc.energy(), pc.position().phi(),
                                 pc.position().z)
                  << endmsg;
          debug() << fmt::format("  - Energy cluster: (E: {}, phi: {}, z: {})", ec.energy(), ec.position().phi(),
                                 ec.position().z)
                  << endmsg;
          debug() << fmt::format("  ---> Merged cluster: (E: {}, phi: {}, z: {})", new_clus.energy(),
                                 new_clus.position().phi(), new_clus.position().z)
                  << endmsg;
        }
      }
    }
    // That's all!

    return StatusCode::SUCCESS;
  }
};
DECLARE_COMPONENT(EnergyPositionClusterMerger)

} // namespace Jug::Reco
