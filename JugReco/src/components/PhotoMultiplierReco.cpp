/*  General PhotoMultiplier Reconstruction
 *
 *  Apply the given quantum efficiency for photon detection
 *  Converts the number of s to signal amplitude
 *
 *  Author: Chao Peng (ANL)
 *  Date: 10/03/2020
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
#include "eicd/PMTHitCollection.h"
#include "eicd/RawPMTHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {
class PhotoMultiplierReco : public GaudiAlgorithm
{
public:
    DataHandle<eic::RawPMTHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::PMTHitCollection>
        m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
    Gaudi::Property<double> m_timeStep{this, "timeStep", 0.0625*ns};
    Gaudi::Property<double> m_minNpe{this, "minNpe", 0.5};
    Gaudi::Property<double> m_speMean{this, "speMean", 80.0};
    Gaudi::Property<double> m_pedMean{this, "pedMean", 200.0};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    PhotoMultiplierReco(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputHitCollection",   m_inputHitCollection,   "");
        declareProperty("outputHitCollection",  m_outputHitCollection,  "");
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
	    const auto &rawhits = *m_inputHitCollection.get();
        // Create output collections
        auto &hits = *m_outputHitCollection.createAndPut();

        // reconstrut number of photo-electrons and time
        for (auto &rh : rawhits) {
            float npe = (rh.amplitude() - m_pedMean)/m_speMean;
            if (npe >= m_minNpe) {
                float time = rh.timeStamp()*m_timeStep;
                auto id = rh.cellID();
                // global positions
                auto gpos = m_geoSvc->cellIDPositionConverter()->position(id);
                // local positions
                auto pos = m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().position();
                hits.push_back(eic::PMTHit{
                    id, npe, time,
                    {gpos.x(), gpos.y(), gpos.z()},
                    {pos.x(), pos.y(), pos.z()}
                });
            }
        }

        return StatusCode::SUCCESS;
    }
};

DECLARE_COMPONENT(PhotoMultiplierReco)

} // namespace Jug::Reco


