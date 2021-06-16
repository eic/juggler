/*
 *  Reconstruct the cluster with Center of Gravity method
 *  Logarithmic weighting is used for mimicing energy deposit in transverse direction
 *
 *  Author: Chao Peng (ANL), 09/27/2020
 */
#include <map>
#include <cstring>
#include <algorithm>
#include <functional>

#include "fmt/format.h"
#include <boost/range/adaptor/map.hpp>
#include "boost/algorithm/string/join.hpp"

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

// weighting functions (with place holders for hit energy, total energy, one parameter and module type enum
static double constWeight(double /*E*/, double /*tE*/, double /*p*/, int /*type*/) { return 1.0; }
static double linearWeight(double E, double /*tE*/, double /*p*/, int /*type*/) { return E; }
static double logWeight(double E, double tE, double base, int /*type*/)
{
    return std::max(0., base + std::log(E/tE));
}

static const std::map<std::string, std::function<double(double, double, double, int)>> weightMethods {
    {"none", constWeight},
    {"linear", linearWeight},
    {"log", logWeight},
};

class ClusterRecoCoG : public GaudiAlgorithm
{
public:
    Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};
    Gaudi::Property<double> m_logWeightBase{this, "logWeightBase", 3.6};
    Gaudi::Property<double> m_depthCorrection{this, "depthCorrection", 0.0};
    Gaudi::Property<std::string> m_energyWeight{this, "energyWeight", "log"};
    Gaudi::Property<std::string> m_moduleDimZName{this, "moduleDimZName", ""};
    DataHandle<eic::ClusterCollection>
        m_clusterCollection{"clusterCollection", Gaudi::DataHandle::Reader, this};
    // Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;
    double m_depthCorr;
    std::function<double(double, double, double, int)> weightFunc;

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

        // select weighting method
        std::string ew = m_energyWeight.value();
        // make it case-insensitive
        std::transform(ew.begin(), ew.end(), ew.begin(), [] (char s) { return std::tolower(s); });
        auto it = weightMethods.find(ew);
        if (it == weightMethods.end()) {
            error() << fmt::format("Cannot find energy weighting method {}, choose one from [{}]",
                                   m_energyWeight,
                                   boost::algorithm::join(weightMethods | boost::adaptors::map_keys, ", "))
                    << endmsg;
            return StatusCode::FAILURE;
        }
        weightFunc = it->second;
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
            float w = weightFunc(hit.energy(), totalE, m_logWeightBase.value(), 0);
            tw += w;
            x += hit.x() * w;
            y += hit.y() * w;
            z += hit.z() * w;
            /*
            debug() << hit.cellID()
                    << "(" << hit.local_x() << ", " << hit.local_y() << ", " << hit.local_z() << "), "
                    << "(" << hit.x() << ", " << hit.y() << ", " << hit.z() << "), "
                    << endmsg;
            */
        }
        res.position({x/tw, y/tw, z/tw});
        // convert global position to local position, use the cell with max edep as a reference
        auto volman = m_geoSvc->detector()->volumeManager();
        auto alignment = volman.lookupDetElement(centerID).nominal();
        auto lpos = alignment.worldToLocal(dd4hep::Position(res.x(), res.y(), res.z()));

        // TODO: may need convert back to have depthCorrection in global positions
        res.local({lpos.x(), lpos.y(), lpos.z() + m_depthCorrection});
        return res;
    }
};

DECLARE_COMPONENT(ClusterRecoCoG)

} // namespace Jug::Reco

