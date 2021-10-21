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
#include "eicd/PMTHitCollection.h"
#include "eicd/RingImageCollection.h" // TODO: remove when no longer needed

using namespace Gaudi::Units;

namespace Jug::Reco {
  /** RICH IRT
   *
   * \ingroup reco
   */
  class TestIRT : public GaudiAlgorithm, AlgorithmIDMixin<> {
  public:

    // local vars and pointers
    SmartIF<IGeoSvc> m_geoSvc; // geometry service
    DataHandle<eic::PMTHitCollection> m_inputHitCollection{
      "inputHitCollection",
      Gaudi::DataHandle::Reader,
      this
    };
    DataHandle<eic::RingImageCollection> m_outputPidCollection{ // TODO: change `RingImageCollection`
      "outputPidCollection",
      Gaudi::DataHandle::Writer,
      this
    };

    // ------------------------------------------------------------

    // constructor
    TestIRT(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc)
      , AlgorithmIDMixin(name, info())
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputPidCollection", m_outputPidCollection, "");
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
      return StatusCode::SUCCESS;
    }


    StatusCode execute() override
    {
      // input collections
      const auto& rawhits = *m_inputHitCollection.get();
      // output collections
      auto& results = *m_outputPidCollection.createAndPut();

      // algorithm
      //auto alg = IRT CODE //// TODO

      return StatusCode::SUCCESS;
    }


    StatusCode finalize() override {
      info() << "TestIRT: Finalizing..." << endmsg;
      return Algorithm::finalize(); // must be executed last
    }


  };

  DECLARE_COMPONENT(TestIRT)

} // namespace Jug::Reco
