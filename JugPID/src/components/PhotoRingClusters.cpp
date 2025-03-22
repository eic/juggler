// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Whitney Armstrong, Sylvester Joosten

/*  Clustering Algorithm for Ring Imaging Cherenkov (RICH) events
 *
 *  Author: Chao Peng (ANL)
 *  Date: 10/04/2020
 *
 */

#include <algorithm>

#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>

// Event Model related classes
#include "FuzzyKClusters.h"
#include "edm4eic/PMTHitCollection.h"
#include "edm4eic/RingImageCollection.h"

using namespace Gaudi::Units;
using namespace Eigen;

namespace Jug::Reco {

/**  Clustering Algorithm for Ring Imaging Cherenkov (RICH) events.
 *
 * \ingroup reco
 */
class PhotoRingClusters : public Gaudi::Algorithm {
private:
  mutable DataHandle<const edm4eic::PMTHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4eic::RingImageCollection> m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer,
                                                                  this};
  // @TODO
  // A more realistic way is to have tracker info as the input to determine how much clusters should be found
  Gaudi::Property<int> m_nRings{this, "nRings", 1};
  Gaudi::Property<int> m_nIters{this, "nIters", 1000};
  Gaudi::Property<double> m_q{this, "q", 2.0};
  Gaudi::Property<double> m_eps{this, "epsilon", 1e-4};
  Gaudi::Property<double> m_minNpe{this, "minNpe", 0.5};
  // Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

public:
  // ill-formed: using Gaudi::Algorithm::GaudiAlgorithm;
  PhotoRingClusters(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputClusterCollection", m_outputClusterCollection, "");
  }

  StatusCode initialize() override {
    if (Gaudi::Algorithm::initialize().isFailure()) {
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

  StatusCode execute(const EventContext&) const override {
    // input collections
    const auto& rawhits = *m_inputHitCollection.get();
    // Create output collections
    auto& clusters = *m_outputClusterCollection.createAndPut();

    // algorithm
    auto alg = fkc::KRings();

    // fill data
    MatrixXd data(rawhits.size(), 2);
    for (int i = 0; i < data.rows(); ++i) {
      if (rawhits[i].getNpe() > m_minNpe) {
        data.row(i) << rawhits[i].getLocal().x, rawhits[i].getLocal().y;
      }
    }

    // clustering
    auto res = alg.Fit(data, m_nRings, m_q, m_eps, m_nIters);

    // local position
    // @TODO: Many fields in RingImage not filled, need to assess
    //        if those are in fact needed
    for (int i = 0; i < res.rows(); ++i) {
      auto cl = clusters.create();
      cl.setPosition({static_cast<float>(res(i, 0)), static_cast<float>(res(i, 1)), 0});
      // @TODO: positionError() not set
      // @TODO: theta() not set
      // @TODO: thetaError() not set
      cl.setRadius(res(i, 2));
      // @TODO: radiusError not set
    }

    return StatusCode::SUCCESS;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(PhotoRingClusters)

} // namespace Jug::Reco
