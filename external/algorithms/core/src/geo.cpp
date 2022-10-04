// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Implementation of the geo service
//
#include <algorithms/geo.h>

namespace algorithms {

void GeoSvc::init(dd4hep::Detector* det) {
  if (det) {
    info() << "Initializing geometry service from pre-initialized detector" << endmsg;
    m_detector = det;
    // no detector given, need to self-initialize
  } else {
    info() << "No external detector provided, self-initializing" << endmsg;
    m_detector = &(dd4hep::Detector::getInstance());
    if (m_xml_list.empty()) {
      // TODO handle error
    }
    for (std::string_view name : m_xml_list) {
      info() << fmt::format("Loading compact file: {}", "name") << endmsg;
      m_detector->fromCompact(std::string(name));
    }
    m_detector->volumeManager();
    m_detector->apply("DD4hepVolumeManager", 0, nullptr);
  }
  // always: instantiate cellIDConverter
  m_converter = std::make_unique<const dd4hep::rec::CellIDPositionConverter>(*m_detector);
}
} // namespace algorithms

