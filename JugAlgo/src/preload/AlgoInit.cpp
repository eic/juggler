// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024, Wouter Deconinck

#include <DD4hep/Detector.h>

#include <algorithms/geo.h>
#include <algorithms/service.h>
#include <algorithms/random.h>
#include <EICrecon/algorithms/interfaces/ParticleSvc.h>

class AlgoInit {
private:
  std::unique_ptr<const dd4hep::Detector> m_detector{nullptr};
public:
  AlgoInit() {
    auto detector = dd4hep::Detector::make_unique("");

    dd4hep::Readout readoutCalorimeter(std::string("MockCalorimeterHits"));
    dd4hep::IDDescriptor id_desc_calorimeter("MockCalorimeterHits", "system:8,layer:8,x:8,y:8");
    dd4hep::Segmentation segmentationCalorimeter("CartesianGridXY","CalorimeterHitsSeg", id_desc_calorimeter.decoder());
    readoutCalorimeter.setIDDescriptor(id_desc_calorimeter);
    readoutCalorimeter.setSegmentation(segmentationCalorimeter);
    detector->add(id_desc_calorimeter);
    detector->add(readoutCalorimeter);

    dd4hep::Readout readoutTracker(std::string("MockTrackerHits"));
    dd4hep::IDDescriptor id_desc_tracker("MockTrackerHits", "system:8,layer:8,x:8,y:8");
    dd4hep::Segmentation segmentationTracker("CartesianGridXY","TrackerHitsSeg", id_desc_tracker.decoder());
    readoutTracker.setIDDescriptor(id_desc_tracker);
    readoutTracker.setSegmentation(segmentationTracker);
    detector->add(id_desc_tracker);
    detector->add(readoutTracker);

    m_detector = std::move(detector);

    auto& serviceSvc = algorithms::ServiceSvc::instance();
    [[maybe_unused]] auto& geoSvc = algorithms::GeoSvc::instance();
    serviceSvc.setInit<algorithms::GeoSvc>([this](auto&& g) {
      g.init(this->m_detector.get());
    });

    [[maybe_unused]] auto& particleSvc = algorithms::ParticleSvc::instance();
    serviceSvc.setInit<algorithms::ParticleSvc>([](auto&& p) {
      p.init(std::make_shared<algorithms::ParticleSvc::ParticleMap>());
    });

    [[maybe_unused]] auto& randomSvc = algorithms::RandomSvc::instance();
    serviceSvc.setInit<algorithms::RandomSvc>([](auto&& r) {
      r.setProperty("seed", static_cast<size_t>(1));
      r.init();
    });

    serviceSvc.init();
  }
};

const AlgoInit algoinit;

