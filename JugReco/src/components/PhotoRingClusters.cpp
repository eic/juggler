/*  Clustering Algorithm for Ring Imaging Cherenkov (RICH) events
 *
 *  Author: Chao Peng (ANL)
 *  Date: 10/04/2020
 *
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
#include "eicd/RIChClusterCollection.h"
#include "FuzzyKClusters.h"

using namespace Gaudi::Units;
using namespace Eigen;


namespace Jug::Reco {
class PhotoRingClusters : public GaudiAlgorithm
{
public:
    DataHandle<eic::PMTHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::RIChClusterCollection>
        m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer, this};
    Gaudi::Property<double> m_minNpe{this, "minNpe", 0.0};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    PhotoRingClusters(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputHitCollection",   m_inputHitCollection,   "");
        declareProperty("outputClusterCollection",  m_outputClusterCollection,  "");
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
        auto &clusters = *m_outputClusterCollection.createAndPut();

        // algorithm
        auto alg = fkc::KRings();

        MatrixXd data(rawhits.size(), 2);
        for (int i = 0; i < data.rows(); ++i) {
            data.row(i) << rawhits[i].local_x(), rawhits[i].local_y();
        }

        return StatusCode::SUCCESS;
    }
};

DECLARE_COMPONENT(PhotoRingClusters)

} // namespace Jug::Reco


