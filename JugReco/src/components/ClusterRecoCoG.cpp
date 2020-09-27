/*
 *  Reconstruct the cluster with Center of Gravity method
 *  Logarithmic weighting is used for mimicing energy deposit in transverse direction
 *
 *  Author: Chao Peng (ANL), 09/27/2020
 */
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
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {
class ClusterRecoCoG : public GaudiAlgorithm
{
public:
    Gaudi::Property<double> m_logWeightBase{this, "logWeightBase", 3.6};
    DataHandle<eic::ClusterCollection>
        m_clusterCollection{"clusterCollection", Gaudi::DataHandle::Reader, this};

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ClusterRecoCoG(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("clusterCollection", m_clusterCollection, "");
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
	    const auto &clusters = *m_clusterCollection.get();
        // reconstruct hit position for the cluster
        for (auto &cl : clusters) {
            reconstruct(cl);
            info() << cl.energy()/GeV << " GeV, (" << cl.position()[0]/mm << ", "
                   << cl.position()[1]/mm << ", " << cl.position()[2]/mm << ")" << endmsg;
        }

        return StatusCode::SUCCESS;
    }

private:
    void reconstruct(eic::Cluster cl) {
        float totalE = 0.;
        for (auto &hit : cl.hits()) {
            totalE += hit.energy();
        }
        cl.energy(totalE);

        // center of gravity with logarithmic weighting
        float totalW = 0., x = 0., y = 0., z = 0.;
        for (auto &hit : cl.hits()) {
            // suppress low energy contributions
            float weight = std::max(0., m_logWeightBase + std::log(hit.energy()/totalE));
            totalW += weight;
            x += hit.position().x * weight;
            y += hit.position().y * weight;
            z += hit.position().z * weight;
        }
        cl.position() = std::array<float, 3>{x/totalW, y/totalW, z/totalW};
    }
};

DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

