/*  General PhotoMultiplier Reconstruction
 *
 *  Estimate the number of photo-electrons and convert timeStamp to time
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

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "eicd/PMTHitCollection.h"
#include "eicd/RawPMTHitCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  /**  General PhotoMultiplier Reconstruction
   *
   *  Estimate the number of photo-electrons and convert timeStamp to time
   *  Collect cell information
   *
   * \ingroup reco
   */
  class PhotoMultiplierReco : public GaudiAlgorithm {
  public:
    DataHandle<eic::RawPMTHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::PMTHitCollection>    m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
    Gaudi::Property<double>              m_timeStep{this, "timeStep", 0.0625 * ns};
    Gaudi::Property<double>              m_minNpe{this, "minNpe", 0.0};
    Gaudi::Property<double>              m_speMean{this, "speMean", 80.0};
    Gaudi::Property<double>              m_pedMean{this, "pedMean", 200.0};
    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    PhotoMultiplierReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                << endmsg;
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& rawhits = *m_inputHitCollection.get();
      // Create output collections
      auto& hits = *m_outputHitCollection.createAndPut();

      // reconstruct number of photo-electrons and time
      for (const auto& rh : rawhits) {
        float npe = (rh.amplitude() - m_pedMean) / m_speMean;
        if (npe >= m_minNpe) {
          float time = rh.time() * (m_timeStep / ns);
          auto  id   = rh.cellID();
          // global positions
          auto gpos = m_geoSvc->cellIDPositionConverter()->position(id);
          // local positions
          auto pos =
              m_geoSvc->cellIDPositionConverter()->findContext(id)->volumePlacement().position();
          // cell dimension
          auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(id);
          hits.push_back(eic::PMTHit{
              rh.cellID(),
              rh.ID(),
              npe,
              time,
              m_timeStep / ns,
              {gpos.x(), gpos.y(), gpos.z()},
              {pos.x(), pos.y(), pos.z()},
              {dim[0]/mm, dim[1]/mm, dim[2]/mm}});
        }
      }

      return StatusCode::SUCCESS;
    }
  };

  DECLARE_COMPONENT(PhotoMultiplierReco)

} // namespace Jug::Reco

