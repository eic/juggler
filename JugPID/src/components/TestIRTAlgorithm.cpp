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

// IRT
#include "IRT/CherenkovDetectorCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {
  /** RICH IRT
   *
   * \ingroup reco
   */
  class TestIRTAlgorithm : public GaudiAlgorithm, AlgorithmIDMixin<> {

  public:

    // local vars and pointers
    SmartIF<IGeoSvc> m_geoSvc; // geometry service
    DataHandle<eic::PMTHitCollection> m_inputHitCollection{
      "inputHitCollection",
      Gaudi::DataHandle::Reader,
      this
    };
    DataHandle<eic::RingImageCollection> m_outputClusterCollection{ // TODO: change `RingImageCollection`
      "outputClusterCollection",
      Gaudi::DataHandle::Writer,
      this
    };
    //Gaudi::Property<int>    m_nRings{this, "nRings", 1}; // TODO: example property

    // ------------------------------------------------------------

    // constructor
    TestIRTAlgorithm(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info())
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputClusterCollection", m_outputClusterCollection, "");
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

      info() << "IRT AQUI-----------------------------------------" << endmsg;
      auto geometry = new CherenkovDetectorCollection();
      geometry->AddNewDetector();
      auto detector = geometry->GetDetector(0);
      info() << "-----------------------------------------" << endmsg;

      return StatusCode::SUCCESS;
    }


    StatusCode execute() override {
      // input collections
      //const auto& rawhits = *m_inputHitCollection.get();
      // output collections
      //auto& results = *m_outputClusterCollection.createAndPut();

      // algorithm
      //auto alg = IRT CODE //// TODO

      return StatusCode::SUCCESS;
    }


    /*
    StatusCode finalize() override {
      info() << "TestIRTAlgorithm: Finalizing..." << endmsg;
      return Algorithm::finalize(); // must be executed last
    }
    */


  };

  DECLARE_COMPONENT(TestIRTAlgorithm)

} // namespace Jug::Reco
