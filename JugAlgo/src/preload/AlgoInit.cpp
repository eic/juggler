// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024, Wouter Deconinck

#include <DD4hep/Detector.h>

#include <algorithms/geo.h>
#include <algorithms/service.h>

class AlgoInit {
private:
  std::unique_ptr<const dd4hep::Detector> m_detector{nullptr};
public:
  AlgoInit() {
    auto detector = dd4hep::Detector::make_unique("");
    dd4hep::Readout readout(std::string("MockCalorimeterHits"));
    dd4hep::IDDescriptor id_desc("MockCalorimeterHits", "system:8,layer:8,x:8,y:8");
    readout.setIDDescriptor(id_desc);
    detector->add(id_desc);
    detector->add(readout);

    dd4hep::Readout readoutTracker(std::string("MockTrackerHits"));
    dd4hep::IDDescriptor id_desc_tracker("MockTrackerHits", "system:8,layer:8,x:8,y:8");
    //Create segmentation with 1x1 mm pixels
    dd4hep::Segmentation segmentation("CartesianGridXY","TrackerHitsSeg", id_desc_tracker.decoder());
    readoutTracker.setIDDescriptor(id_desc_tracker);
    readoutTracker.setSegmentation(segmentation);
    detector->add(id_desc_tracker);
    detector->add(readoutTracker);

    m_detector = std::move(detector);

    auto& serviceSvc = algorithms::ServiceSvc::instance();
    [[maybe_unused]] auto& geoSvc = algorithms::GeoSvc::instance();
    serviceSvc.setInit<algorithms::GeoSvc>([this](auto&& g) {
      g.init(this->m_detector.get());
    });

    serviceSvc.init();
  }
};

const AlgoInit algoinit;

