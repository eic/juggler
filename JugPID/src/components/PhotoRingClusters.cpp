/*  Clustering Algorithm for Ring Imaging Cherenkov (RICH) events
 *
 *  Author: Chao Peng (ANL)
 *  Date: 10/04/2020
 *
 */

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
#include "FuzzyKClusters.h"
#include "eicd/PMTHitCollection.h"
#include "eicd/RingImageCollection.h"

using namespace Gaudi::Units;
using namespace Eigen;

namespace Jug::Reco {

  /**  Clustering Algorithm for Ring Imaging Cherenkov (RICH) events.
   *
   * \ingroup reco
   */
  class PhotoRingClusters : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:
    DataHandle<eic::PMTHitCollection>      m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::RingImageCollection> m_outputClusterCollection{"outputClusterCollection",
                                                                     Gaudi::DataHandle::Writer, this};
    // @TODO
    // A more realistic way is to have tracker info as the input to determine how much clusters should be found
    Gaudi::Property<int>    m_nRings{this, "nRings", 1};
    Gaudi::Property<int>    m_nIters{this, "nIters", 1000};
    Gaudi::Property<double> m_q{this, "q", 2.0};
    Gaudi::Property<double> m_eps{this, "epsilon", 1e-4};
    Gaudi::Property<double> m_minNpe{this, "minNpe", 0.5};
    // Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    PhotoRingClusters(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin(name, info())
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputClusterCollection", m_outputClusterCollection, "");
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
      const auto& rawhits = *m_inputHitCollection.get();
      // Create output collections
      auto& clusters = *m_outputClusterCollection.createAndPut();

      // algorithm
      auto alg = fkc::KRings();

      // fill data
      MatrixXd data(rawhits.size(), 2);
      for (int i = 0; i < data.rows(); ++i) {
        if (rawhits[i].npe() > m_minNpe) {
          data.row(i) << rawhits[i].local().x, rawhits[i].local().y;
        }
      }

      // clustering
      auto res = alg.Fit(data, m_nRings, m_q, m_eps, m_nIters);

      // local position
      // @TODO: Many fields in RingImage not filled, need to assess
      //        if those are in fact needed
      for (int i = 0; i < res.rows(); ++i) {
        auto cl = clusters.create();
        cl.ID({i, algorithmID()});
        cl.position({res(i, 0), res(i, 1), 0});
        // @TODO: positionError() not set
        // @TODO: theta() not set
        // @TODO: thetaError() not set
        cl.radius(res(i, 2));
        // @TODO: radiusError not set
      }

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(PhotoRingClusters)

} // namespace Jug::Reco
