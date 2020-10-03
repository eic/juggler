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
    Gaudi::Property<std::string> m_moduleDimZName{this, "moduleDimZName", "CrystalBox_z_length"};
    DataHandle<eic::ClusterCollection>
        m_clusterCollection{"clusterCollection", Gaudi::DataHandle::Reader, this};
    // Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;
    double m_depthCorr;

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
        m_geoSvc = service("GeoSvc");
        if (!m_geoSvc) {
            error() << "Unable to locate Geometry Service. "
                    << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
            return StatusCode::FAILURE;
        }
	// depth: z length of the crystal block
        m_depthCorr = m_geoSvc->detector()->constantAsDouble(m_moduleDimZName);
        //info() << "z_length " << depth << endmsg;
        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // input collections
        auto &clusters = *m_clusterCollection.get();
        // reconstruct hit position for the cluster
        for (auto &cl : clusters) {
            reconstruct(cl);
            // info() << cl.energy()/GeV << " GeV, (" << cl.position().x/mm << ", "
            //        << cl.position().y/mm << ", " << cl.position().z/mm << ")" << endmsg;
        }

        return StatusCode::SUCCESS;
    }

private:
    void reconstruct(eic::Cluster cl)
    {
        // no hits
        if (cl.hits_size() == 0) {
            return;
        }

        // calculate total energy, find the cell with the maximum energy deposit
        float totalE = 0., maxE = 0.;
        auto centerID = cl.hits_begin()->cellID();
        for (auto &hit : cl.hits()) {
            auto energy = hit.energy();
            totalE += energy;
            if (energy > maxE) {
                maxE = energy;
                centerID = hit.cellID();
            }
        }
        cl.energy(totalE);

        // center of gravity with logarithmic weighting
        float tw = 0., x = 0., y = 0., z = 0.;
        for (auto &hit : cl.hits()) {
            // suppress low energy contributions
            float w = std::max(0., m_logWeightBase + std::log(hit.energy()/totalE));
            tw += w;
            x += hit.local_x() * w;
            y += hit.local_y() * w;
            z += hit.local_z() * w;
        }

        // convert local position to global position, use the cell with max edep as a reference
        auto volman = m_geoSvc->detector()->volumeManager();
        auto alignment = volman.lookupDetector(centerID).nominal();
        auto gpos = alignment.localToWorld(dd4hep::Position(x/tw, y/tw, z/tw + m_depthCorr));

        cl.position({gpos.x(), gpos.y(), gpos.z()});
    }
};

DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

