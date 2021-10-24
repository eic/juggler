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

#include "eicd/PMTHitCollection.h"
#include "eicd/RingImageCollection.h" // TODO: remove when no longer needed

#include "IRT/ParametricSurface.h"
#include "IRT/CherenkovRadiator.h"
#include "IRT/OpticalBoundary.h"
#include "IRT/CherenkovDetectorCollection.h"
#include "IRT/CherenkovPhotonDetector.h"

#include "TVector3.h"

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

    // constructor
    TestIRTAlgorithm(const std::string& name, ISvcLocator* svcLoc) 
      : GaudiAlgorithm(name, svcLoc), AlgorithmIDMixin(name, info())
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputClusterCollection", m_outputClusterCollection, "");
    }


    // INITIALIZE .....................................................
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

      // create IRT detector collection geometry
      auto irtGeo = new CherenkovDetectorCollection();
      irtGeo->AddNewDetector();
      auto irtDet = irtGeo->GetDetector(0);
      // TODO: copy over Alexander's FIXME comments

      // get detector elements
      auto dd4Det = m_geoSvc->detector();
      auto erichDet = dd4Det->detector("ERICH");
      //auto drichDet = dd4Det->detector("DRICH");

      // get detector positions
      auto erichPos = erichDet.placement().position();
      //auto drichPos = drichDet.placement().position();

      // set IRT container volume
      auto boundary = new FlatSurface(TVector3(0,0,0), TVector3(1,0,0), TVector3(0,1,0));
      irtGeo->SetContainerVolume(irtDet, (G4LogicalVolume*)(0x0), 0, boundary);

      // eRICH loop :::::::::::::::::::::::
      int sector;
      for(auto const& [erichDEname, erichDE] : erichDet.children()) {
        //info() << "FOUND ERICH DE: " << erichDEname << endmsg;
        auto pos = erichPos + erichDE.placement().position();

        // aerogel and filter
        if(erichDEname.find("aerogel")!=std::string::npos || erichDEname.find("filter")!=std::string::npos) {
          double thickness = 2 * erichDE.volume().boundingBox().dimensions()[2];
          sector = erichDE.id();
          info() << "ERICH RADIATOR: " << erichDEname
            << "\n\t(x,y,z)-position = " << pos.x() << ", " << pos.y() << ", " << pos.z()
            << "\n\tsector = " << sector
            << "\n\tthickness = " << thickness
            << endmsg;
          if(sector==0) {
            if(erichDEname.find("aerogel")!=std::string::npos) {
              auto aerogelSurf = new FlatSurface(
                  (1/mm)*TVector3(0,0,pos.z()), // TODO: correct position?
                  TVector3(1,0,0),
                  TVector3(0,1,0)
                  );
              irtGeo->AddFlatRadiator(irtDet, (G4LogicalVolume*)(0x1), 0, aerogelSurf, thickness/mm); // TODO: correct units?
            } else {
              auto filterSurf = new FlatSurface(
                  (1/mm)*TVector3(0,0,pos.z()-0.01), // TODO: correct position?
                  TVector3(1,0,0),
                  TVector3(0,1,0)
                  );
              irtGeo->AddFlatRadiator(irtDet, (G4LogicalVolume*)(0x2), 0, filterSurf, thickness/mm);
            }
          }
        } // end if aerogel


        // sensors
        if(erichDEname.find("sensor")!=std::string::npos) {
          info() << "ERICH SENSOR: " << erichDEname
            << "\n\t(x,y,z)-position = " << pos.x() << ", " << pos.y() << ", " << pos.z()
            << "\n\tid = " << erichDE.id()
            << endmsg;
          auto sensorSurf = new FlatSurface(
              (1/mm)*TVector3(0.0, 0.0, pos.z()), // TODO: why not `pos.x(), pos.y(), pos.z()`?
              TVector3(1,0,0),
              TVector3(0,1,0)
              );
          irtDet->AddPhotonDetector(new CherenkovPhotonDetector(0,0,sensorSurf));
        }

      } // end loop over eRICH detector elements


      // set radiator refractive indices
      info() << "C4F10 RINDEX: " << endmsg;
      auto gasRmatrix = dd4Det->material("C4F10_ERICH").property("RINDEX"); // TODO: get material from gas volume instead
      for(size_t r=0; r<gasRmatrix->GetRows(); r++) {
        info() << "   " << gasRmatrix->Get(r,0) << "   " << gasRmatrix->Get(r,1) << endmsg;
      };
      // TODO: irtDet->Radiators()[0]->SetReferenceRefractiveIndex( N ) // aerogel
      // TODO: irtDet->Radiators()[1]->SetReferenceRefractiveIndex( N ) // filter
      // TODO: irtDet->Radiators()[2]->SetReferenceRefractiveIndex( N ) // gas

      return StatusCode::SUCCESS;
    }


    // EXECUTE .....................................................
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
