#include <cmath>
// Gaudi
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/ToolHandle.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugTrack/ProtoTrack.hpp"
#include "JugTrack/Track.hpp"

#include "Math/Vector3D.h"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

/** Template finder.
 *
 *  \ingroup tracking
 */
class FinderAlgoTemplate : public GaudiAlgorithm {
public:
  DataHandle<eicd::TrackerHitCollection> m_inputTrackerHits{"inputTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<Jug::ProtoTrackContainer> m_outputProtoTracks{"outputProtoTracks", Gaudi::DataHandle::Writer, this};

public:
  FinderAlgoTemplate(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputTrackerHits", m_inputTrackerHits, "tracker hits whose indices are used in proto-tracks");
    declareProperty("outputProtoTracks", m_outputProtoTracks, "grouped hit indicies");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collection
    //const eicd::TrackerHitCollection* hits = m_inputTrackerHits.get();
    // Create output collections
    /*auto proto_tracks = */m_outputProtoTracks.createAndPut();

    return StatusCode::SUCCESS;
  }
};
DECLARE_COMPONENT(FinderAlgoTemplate)

} // namespace Jug::Reco
