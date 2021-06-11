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
    Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.6};
    Gaudi::Property<double> m_logWeightBase{this, "logWeightBase", 3.6};
    Gaudi::Property<double> m_depthCorrection{this, "depthCorrection", 0.0};
    Gaudi::Property<std::string> m_moduleDimZName{this, "moduleDimZName", ""};
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
        // update depth correction if a name is provided
        if (!m_moduleDimZName.value().empty()) {
            m_depthCorrection = m_geoSvc->detector()->constantAsDouble(m_moduleDimZName);
        }
        //info() << "z_length " << depth << endmsg;
        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // input collections
        auto &clusters = *m_clusterCollection.get();
        // reconstruct hit position for the cluster
        for (auto cl : clusters) {
            auto hit = reconstruct(cl);
            cl.nhits(cl.hits_size());
            cl.edep(hit.energy());
            cl.energy(hit.energy()/m_sampFrac);
            cl.position(hit.position());
            cl.polar(cart_to_polar(hit.position()));
            debug() << cl.hits_size() << " hits: " << cl.energy()/GeV << " GeV, (" << cl.position().x/mm << ", "
                    << cl.position().y/mm << ", " << cl.position().z/mm << ")" << endmsg;
        }

        return StatusCode::SUCCESS;
    }

private:
    template<typename T1>
    eic::VectorPolar cart_to_polar(const T1 &cart) {
        auto r = std::sqrt(cart.x*cart.x + cart.y*cart.y + cart.z*cart.z);
        return eic::VectorPolar{r, std::acos(cart.z/r), std::atan2(cart.y, cart.x)};
    }

    eic::CalorimeterHit reconstruct(eic::Cluster cl) const
    {
        eic::CalorimeterHit res;
        // no hits
        if (cl.hits_size() == 0) {
            return res;;
        }

        // calculate total energy, find the cell with the maximum energy deposit
        float totalE = 0., maxE = 0.;
        auto centerID = cl.hits_begin()->cellID();
        for (auto &hit : cl.hits()) {
            // info() << "hit energy = " << hit.energy() << endmsg;
            auto energy = hit.energy();
            totalE += energy;
            if (energy > maxE) {
                maxE = energy;
                centerID = hit.cellID();
            }
        }
        res.cellID(centerID);
        res.energy(totalE);

        // center of gravity with logarithmic weighting
        float tw = 0., x = 0., y = 0., z = 0.;
        for (auto &hit : cl.hits()) {
            // suppress low energy contributions
            // info() << std::log(hit.energy()/totalE) << endmsg;
            float w = std::max(0., m_logWeightBase + std::log(hit.energy()/totalE));
            tw += w;
            x += hit.local_x() * w;
            y += hit.local_y() * w;
            z += hit.local_z() * w;
            /*
            debug() << hit.cellID()
                    << "(" << hit.local_x() << ", " << hit.local_y() << ", " << hit.local_z() << "), "
                    << "(" << hit.x() << ", " << hit.y() << ", " << hit.z() << "), "
                    << endmsg;
            */
        }
        res.local({x/tw, y/tw, z/tw + m_depthCorrection});

        // convert local position to global position, use the cell with max edep as a reference
        auto volman = m_geoSvc->detector()->volumeManager();
        auto alignment = volman.lookupDetElement(centerID).nominal();
        auto gpos = alignment.localToWorld(dd4hep::Position(res.local_x(), res.local_y(), res.local_z()));

        res.position({gpos.x(), gpos.y(), gpos.z()});
        return res;
    }
};

DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

