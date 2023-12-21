// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Chao Peng, Sylvester Joosten, Whitney Armstrong, Chao Peng

/*  General PhotoMultiplier Reconstruction
 *
 *  Estimate the number of photo-electrons and convert getTimeStamp to time
 *  Collect cell information
 *
 *  Author: Chao Peng (ANL)
 *  Date: 10/03/2020
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

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>

// Event Model related classes
#include "edm4eic/PMTHitCollection.h"
#include "edm4eic/RawTrackerHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

/**  General PhotoMultiplier Reconstruction
 *
 *  Estimate the number of photo-electrons and convert getTimeStamp to time
 *  Collect cell information
 *
 * \ingroup reco
 */
class PhotoMultiplierReco : public GaudiAlgorithm {
private:
  DataHandle<edm4eic::RawTrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::PMTHitCollection> m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
  Gaudi::Property<double> m_timeStep{this, "timeStep", 0.0625 * ns};
  Gaudi::Property<double> m_minNpe{this, "minNpe", 0.0};
  Gaudi::Property<double> m_speMean{this, "speMean", 80.0};
  Gaudi::Property<double> m_pedMean{this, "pedMean", 200.0};
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  std::shared_ptr<const dd4hep::rec::CellIDPositionConverter> m_converter;

public:
  // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
  PhotoMultiplierReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputHitCollection", m_inputHitCollection, "");
    declareProperty("outputHitCollection", m_outputHitCollection, "");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    m_converter = std::make_shared<const dd4hep::rec::CellIDPositionConverter>(*(m_geoSvc->getDetector()));
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& rawhits = *m_inputHitCollection.get();
    // Create output collections
    auto& hits = *m_outputHitCollection.createAndPut();

    // reconstruct number of photo-electrons and time
    for (const auto& rh : rawhits) {
      float npe = (rh.getCharge() - m_pedMean) / m_speMean;
      if (npe >= m_minNpe) {
        float time = rh.getTimeStamp() * (static_cast<float>(m_timeStep) / ns);
        auto id    = rh.getCellID();
        // global positions
        auto gpos = m_converter->position(id);
        // local positions
        auto pos = m_converter->findContext(id)->volumePlacement().position();
        // cell dimension
        auto dim = m_converter->cellDimensions(id);
        hits.push_back(edm4eic::PMTHit{
            rh.getCellID(),
            npe,
            time,
            static_cast<float>(m_timeStep / ns),
            {static_cast<float>(gpos.x()), static_cast<float>(gpos.y()), static_cast<float>(gpos.z())},
            {static_cast<float>(dim[0] / mm), static_cast<float>(dim[1] / mm), static_cast<float>(dim[2] / mm)},
            0, // @FIXME: Add sector similar to CalorimeterHit
            {static_cast<float>(pos.x()), static_cast<float>(pos.y()), static_cast<float>(pos.z())}});
      }
    }

    return StatusCode::SUCCESS;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(PhotoMultiplierReco)

} // namespace Jug::Reco
