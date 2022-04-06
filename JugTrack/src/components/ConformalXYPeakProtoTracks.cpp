// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten

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

#include "TH1F.h"
#include "Math/Vector3D.h"
#include "Math/Vector2D.h"

#include "eicd/TrackerHitCollection.h"

namespace Jug::Reco {

/** Conformal XY hits.
 *
 *  Conformal mapping which turns circles to lines.
 *
 *  \ingroup tracking
 */
class ConformalXYPeakProtoTracks : public GaudiAlgorithm {
private:
  DataHandle<eicd::TrackerHitCollection> m_inputTrackerHits{"inputTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<Jug::ProtoTrackContainer> m_outputProtoTracks{"outputProtoTracks", Gaudi::DataHandle::Writer, this};
  DataHandle<int> m_nProtoTracks{"nProtoTracks", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<int> m_nPhiBins{this, "nPhiBins", 100};

  using ConformalHit = ROOT::Math::XYVector;

public:
  ConformalXYPeakProtoTracks(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputTrackerHits", m_inputTrackerHits, "tracker hits whose indices are used in proto-tracks");
    declareProperty("outputProtoTracks", m_outputProtoTracks, "grouped hit indicies");
    declareProperty("nProtoTracks", m_nProtoTracks, "number of proto tracks");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collection
    const eicd::TrackerHitCollection* hits = m_inputTrackerHits.get();
    // Create output collections
    auto* proto_tracks = m_outputProtoTracks.createAndPut();
    int n_proto_tracks = 0;

    std::vector<ConformalHit> conformal_hits;

    ConformalHit ref_hit(0.0,0.0); // future versions will improve on this.

    TH1F h_phi("h_phi",";phi",m_nPhiBins.value(),-M_PI,M_PI);

    // 1. conformal XY transform hits
    // 2. fill histogram with phi
    for(const auto& ahit : *hits) {
      double xc = ahit.getPosition().x - ref_hit.x();
      double yc = ahit.getPosition().y - ref_hit.y();
      double r = std::hypot(xc, yc);
      conformal_hits.emplace_back(2.0*xc/r,2.0*yc/r);
      h_phi.Fill(conformal_hits.back().phi());
    }
    // 3. Get location of maxima
    std::vector<int> max_bins;
    while(max_bins.size() < 100) {
      int    max_bin = h_phi.GetMaximumBin();
      double max_val = h_phi.GetMaximum();
      if(max_val < 3)  {
        break;
      }
      max_bins.push_back(max_bin);
      h_phi.SetBinContent(max_bin, 0.0); // zero bin and continue
    }
    n_proto_tracks = max_bins.size();
    if (msgLevel(MSG::DEBUG)) {
      debug() << " Found " << n_proto_tracks << " proto tracks." << endmsg;
    }
    // 4. Group hits peaked in phi
    for(auto b : max_bins) {
      Jug::ProtoTrack proto_track; // this is just a std::vector<int>
      for(size_t ihit = 0 ; ihit< hits->size() ; ihit++) {
        double phi = conformal_hits[ihit].phi();
        double bin_phi = h_phi.GetXaxis()->GetBinCenter(b);
        double bin_width = h_phi.GetXaxis()->GetBinWidth(b); /// \todo make bin width an algo parameter
        if (std::abs(phi - bin_phi) < bin_width/2.0) {
          proto_track.push_back(ihit);
        }
      }
      if(proto_track.size()>3) {
        proto_tracks->push_back(proto_track);
      }
    }
    m_nProtoTracks.put(new int{n_proto_tracks});
    // 5. profit!

    return StatusCode::SUCCESS;
  }
};
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(ConformalXYPeakProtoTracks)

} // namespace Jug::Reco
