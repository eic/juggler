// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#include <cmath>
// Gaudi
#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include "Math/Vector3D.h"

#include "edm4eic/TrackerHitCollection.h"

namespace Jug::Reco {

/** Template finder.
 *
 *  \ingroup tracking
 */
class ProtoTrackMatching : public Gaudi::Algorithm {
private:
  mutable DataHandle<edm4eic::TrackerHitCollection> m_inputTrackerHits{"inputTrackerHits", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<ActsExamples::TrackParametersContainer> m_initialTrackParameters{"initialTrackParameters", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<ActsExamples::ProtoTrackContainer> m_inputProtoTracks{"inputProtoTracks", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<ActsExamples::ProtoTrackContainer> m_outputProtoTracks{"matchedProtoTracks", Gaudi::DataHandle::Writer, this};

public:
  ProtoTrackMatching(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
    declareProperty("inputTrackerHits",     m_inputTrackerHits,     "");
    declareProperty("initialTrackParameters", m_initialTrackParameters, "");
    declareProperty("inputProtoTracks",       m_inputProtoTracks,       "");
    declareProperty("matchedProtoTracks",     m_outputProtoTracks, "proto tracks matched to initial track parameters");
  }

  StatusCode initialize() override {
    if (Gaudi::Algorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    // input collection
    
    // hits is unused, commented out for now to avoid compiler warning
    //const auto* const hits              = m_inputTrackerHits.get();
    
    const auto* const proto_tracks      = m_inputProtoTracks.get();
    const auto* const initialParameters = m_initialTrackParameters.get();

    // Create output collections
    auto* matched_proto_tracks = m_outputProtoTracks.createAndPut();

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
      // track_param and iproto_best currently unused, but will be in the future
      // commented out to remove compiler warning
      //const auto& track_param = (*initialParameters)[itrack];
      //int iproto_best = 0;

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
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ProtoTrackMatching)

} // namespace Jug::Reco
