#ifndef JUGGLER_JUGRECO_GenFitTrackFitter_HH
#define JUGGLER_JUGRECO_GenFitTrackFitter_HH

#include "JugTrack/GeometryContainers.hpp"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

#include <functional>
#include <stdexcept>
#include <vector>

#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/BField.h"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/Trajectories.hpp"
#include "JugTrack/ProtoTrack.hpp"

#include "eicd/TrackerHitCollection.h"

#include <random>
#include <stdexcept>

namespace Jug::Reco {

  /** Genfit based tracking algorithm.
   *
   * \ingroup track
   * \ingroup tracker
   */
  class GenFitTrackFitter : public GaudiAlgorithm {
  public:

  public:
    DataHandle<eic::TrackerHitCollection>    m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<TrackParametersContainer> m_initialTrackParameters{"initialTrackParameters", Gaudi::DataHandle::Reader, this};
    DataHandle<ProtoTrackContainer>      m_inputProtoTracks{"inputProtoTracks", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoriesContainer>    m_foundTracks{"foundTracks", Gaudi::DataHandle::Reader, this};
    DataHandle<TrajectoriesContainer>    m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};

    SmartIF<IGeoSvc>                      m_geoSvc;
    //Acts::GeometryContext                 m_geoctx;
    //Acts::CalibrationContext              m_calibctx;
    //Acts::MagneticFieldContext            m_fieldctx;


    GenFitTrackFitter(const std::string& name, ISvcLocator* svcLoc);

    StatusCode initialize() override;
    StatusCode execute() override;
  };


} // namespace Jug::Reco

#endif
