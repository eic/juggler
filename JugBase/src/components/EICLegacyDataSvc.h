// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Benedikt Hegner

#include "EICLegacyDataSvc.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
DECLARE_COMPONENT(EICLegacyDataSvc)

/// Standard Constructor
EICLegacyDataSvc::EICLegacyDataSvc(const std::string& name, ISvcLocator* svc) : PodioLegacyDataSvc(name, svc) {
  declareProperty("inputs", m_filenames = {}, "Names of the files to read");
  declareProperty("input", m_filename = "", "Name of the file to read");
  declareProperty("FirstEventEntry", m_1stEvtEntry = 0, "First event to read");
}

/// Standard Destructor
EICLegacyDataSvc::~EICLegacyDataSvc() {}
