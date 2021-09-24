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
class ProtoTrackMatching : public GaudiAlgorithm {
public:
  DataHandle<eic::TrackerHitCollection> m_inputTrackerHits{"inputTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<TrackParametersContainer> m_initialTrackParameters{"initialTrackParameters", Gaudi::DataHandle::Reader, this};
  DataHandle<ProtoTrackContainer> m_inputProtoTracks{"inputProtoTracks", Gaudi::DataHandle::Reader, this};
  DataHandle<ProtoTrackContainer> m_outputProtoTracks{"matchedProtoTracks", Gaudi::DataHandle::Writer, this};

public:
  ProtoTrackMatching(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputTrackerHits",     m_inputTrackerHits,     "");
    declareProperty("initialTrackParameters", m_initialTrackParameters, "");
    declareProperty("inputProtoTracks",       m_inputProtoTracks,       "");
    declareProperty("matchedProtoTracks",     m_outputProtoTracks, "proto tracks matched to initial track parameters");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collection
    const eic::TrackerHitCollection* hits              = m_inputTrackerHits.get();
    const ProtoTrackContainer*       proto_tracks      = m_inputProtoTracks.get();
    const TrackParametersContainer*  initialParameters = m_initialTrackParameters.get();

    // Create output collections
    auto matched_proto_tracks = m_outputProtoTracks.createAndPut();

    int n_tracks = initialParameters->size();
    int n_proto_tracks = proto_tracks->size();

    // Assuming init track parameters have been match with proto tracks by index
    if(n_proto_tracks <  n_tracks) {
      warning() << " Number of proto tracks does not match the initial track parameters." << endmsg;
      if(n_proto_tracks == 0 ) {
        warning() << " Zero proto tracks to match! " << endmsg;
        return StatusCode::SUCCESS;
      }
    }


    if (msgLevel(MSG::DEBUG)) {
      debug() << " ntracks        " << n_tracks <<  endmsg;
      debug() << " n_proto_tracks " << n_proto_tracks <<  endmsg;
    }
    // do better matching
    for (int itrack = 0; itrack < n_tracks; itrack++) {
      const auto& track_param = (*initialParameters)[itrack];
      int iproto_best = 0;

      /// \todo FIXME
      /// find the best matching proto track:
      /// here we just take the first.
      if( n_proto_tracks <= itrack) {
          matched_proto_tracks->push_back( proto_tracks->back() );
      } else {
          matched_proto_tracks->push_back( (*proto_tracks)[itrack] );
      }
      //for(const auto& aproto : (*protoTracks)){
      //  const auto& proto_track = (*protoTracks)[iproto];
      //  if(n_proto_tracks <  n_tracks) {
      //  }
      //}
    }
    return StatusCode::SUCCESS;
  }
};
DECLARE_COMPONENT(ProtoTrackMatching)

} // namespace Jug::Reco
