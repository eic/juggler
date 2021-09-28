#include <algorithm>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
//#include "GaudiKernel/ToolHandle.h"
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

#include "DD4hep/DD4hepUnits.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
//#include "GaudiExamples/MyTrack.h"
#include "eicd/RawTrackerHitCollection.h"
#include "eicd/TrackerHitCollection.h"

#include "DD4hep/DD4hepUnits.h"

namespace {
inline double get_resolution(const double pixel_size) {
  constexpr const double sqrt_12 = 3.4641016151;
  return pixel_size / sqrt_12;
}
inline double get_variance(const double pixel_size) {
  const double res = get_resolution(pixel_size);
  return res * res;
}
} // namespace

namespace Jug {
  namespace Reco {

    /** Tracker hit reconstruction.
     *
     * \ingroup reco
     */
    class TrackerHitReconstruction : public GaudiAlgorithm, AlgorithmIDMixin<> {
    public:
      Gaudi::Property<double>                  m_timeResolution{this, "timeResolution", 10};
      Rndm::Numbers                            m_gaussDist;
      DataHandle<eic::RawTrackerHitCollection> m_inputHitCollection{
          "inputHitCollection", Gaudi::DataHandle::Reader, this};
      DataHandle<eic::TrackerHitCollection> m_outputHitCollection{"outputHitCollection",
                                                                  Gaudi::DataHandle::Writer, this};
      /// Pointer to the geometry service
      SmartIF<IGeoSvc> m_geoSvc;

    public:
      //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
      TrackerHitReconstruction(const std::string& name, ISvcLocator* svcLoc)
          : GaudiAlgorithm(name, svcLoc)
          , AlgorithmIDMixin<>(name, info())
      {
        declareProperty("inputHitCollection", m_inputHitCollection, "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
      }

      StatusCode initialize() override
      {
        if (GaudiAlgorithm::initialize().isFailure())
          return StatusCode::FAILURE;
        m_geoSvc = service("GeoSvc");
        if (!m_geoSvc) {
          error() << "Unable to locate Geometry Service. "
                  << "Make sure you have GeoSvc and SimSvc in the right order in the configuration."
                  << endmsg;
          return StatusCode::FAILURE;
        }
        IRndmGenSvc* randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        StatusCode sc = m_gaussDist.initialize(randSvc, Rndm::Gauss(0.0, m_timeResolution.value()));
        if (!sc.isSuccess()) {
          return StatusCode::FAILURE;
        }
        return StatusCode::SUCCESS;
      }

      StatusCode execute() override
      {
        constexpr auto mm = dd4hep::mm;
        // input collection
        const eic::RawTrackerHitCollection* rawhits = m_inputHitCollection.get();
        // Create output collections
        auto rec_hits = m_outputHitCollection.createAndPut();

        debug() << " raw hits size : " <<  std::size(*rawhits) << endmsg;
        for (const auto& ahit : *rawhits) {
          // debug() << "cell ID : " << ahit.cellID() << endmsg;
          auto pos = m_geoSvc->cellIDPositionConverter()->position(ahit.cellID());
          auto dim = m_geoSvc->cellIDPositionConverter()->cellDimensions(ahit.cellID());

          if (msgLevel(MSG::VERBOSE)) {
            size_t i = 0;
            for (const auto& p : {pos.x(), pos.y(), pos.z()}) {
              verbose() << "position " << i++ << " [mm]: " << p / mm << endmsg;
            }
            verbose() << "dimension size: " << std::size(dim) << endmsg;
            for (size_t j = 0; j < std::size(dim); ++j) {
              verbose() << " - dimension " << j << " size: " << dim[j] << endmsg;
            }
          }
          // Note about variance:
          //    The variance is used to obtain a diagonal covariance matrix.
          //    Note that the covariance matrix is written in DD4hep surface coordinates,
          //    *NOT* global position coordinates. This implies that:
          //      - XY segmentation: xx -> sigma_x, yy-> sigma_y, zz -> 0, tt -> 0
          //      - XZ segmentation: xx -> sigma_x, yy-> sigma_z, zz -> 0, tt -> 0
          //      - XYZ segmentation: xx -> sigma_x, yy-> sigma_y, zz -> sigma_z, tt -> 0
          //    This is properly in line with how we get the local coordinates for the hit
          //    in the TrackerSourceLinker.
          eic::TrackerHit hit{{ahit.ID().value, algorithmID()}, // Hit ID
                              ahit.cellID(),                    // Raw DD4hep cell ID
                              {pos.x() / mm, pos.y() / mm, pos.z() / mm, (float)ahit.time() / 1000}, // mm, ns
                              {get_variance(dim[0] / mm), get_variance(dim[1] / mm), // variance (see note above)
                               std::size(dim) > 2 ? get_variance(dim[2] / mm) : 0., 0},
                              static_cast<float>(ahit.charge() / 1.0e6), // Collected energy (GeV)
                              0.0f};                                     // Error on the energy
          rec_hits->push_back(hit);
        }
        return StatusCode::SUCCESS;
      }
    };
    DECLARE_COMPONENT(TrackerHitReconstruction)

  } // namespace Reco
} // namespace Jug

