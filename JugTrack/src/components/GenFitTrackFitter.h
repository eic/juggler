// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#ifndef JUGGLER_JUGRECO_GenFitTrackFitter_HH
#define JUGGLER_JUGRECO_GenFitTrackFitter_HH

#include <functional>
#include <stdexcept>
#include <vector>

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/BField/DD4hepBField.h"
#include "JugTrack/GeometryContainers.hpp"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Track.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/Trajectories.hpp"
#include "JugTrack/ProtoTrack.hpp"

#include "edm4eic/TrackerHitCollection.h"
#include "edm4eic/TrajectoryCollection.h"
#include "edm4eic/TrackParametersCollection.h"

//genfitk
#include "FieldManager.h"

#include <random>
#include <stdexcept>

namespace Jug::Reco {


  /** Genfit based tracking algorithm.
   *
   * \ingroup tracking
   */
  class GenFitTrackFitter : public GaudiAlgorithm {
  public:

  class FieldImp : public genfit::AbsBField {
  protected:
    dd4hep::Detector* m_detector;
  public:
    FieldImp(dd4hep::Detector* det): m_detector(det) {}
    virtual ~FieldImp() {}

    /** Get the magneticField [kGauss] at position.
     *
     *  Note that tgeo units are used. [cm] and [kGauss].
     */
    TVector3 get(const TVector3& position) const override {
      double pos[3] = {position.x(), position.y(), position.z()};
      double field[3];
      this->get(pos[0], pos[1], pos[2], field[0], field[1], field[2]);
      return {field[0], field[1], field[2]};
    }

    /** Get the magneticField [kGauss] at position.
     *
     *  Note that tgeo units are used. [cm] and [kGauss].
     */
    void get(const double& posX, const double& posY, const double& posZ,
             double& Bx, double& By, double& Bz) const override {
    dd4hep::Position pos(posX,posY,posZ);
    auto field = m_detector->field().magneticField(pos) * (dd4hep::kilogauss / dd4hep::tesla);
    Bx = field.x();
    By = field.y();
    Bz = field.z();
    //return {field.x(), field.y(),field.z()};
    }
  };

public:
  DataHandle<edm4eic::TrackerHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<TrackParametersContainer>  m_initialTrackParameters{"initialTrackParameters", Gaudi::DataHandle::Reader, this};
  DataHandle<ProtoTrackContainer>       m_inputProtoTracks{"inputProtoTracks", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::TrackParametersCollection> m_foundTracks{"trackParameters", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4eic::TrajectoryCollection> m_outputTrajectories{"outputTrajectories", Gaudi::DataHandle::Writer, this};

  SmartIF<IGeoSvc> m_geoSvc;
  // Acts::GeometryContext                 m_geoctx;
  // Acts::CalibrationContext              m_calibctx;
  // Acts::MagneticFieldContext            m_fieldctx;

  std::map<int64_t, std::shared_ptr<genfit::DetPlane>> m_detPlaneMap;
  std::map<int64_t, dd4hep::rec::Surface*> m_surfaceMap;

  GenFitTrackFitter(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode execute() override;
  };


} // namespace Jug::Reco

#endif
