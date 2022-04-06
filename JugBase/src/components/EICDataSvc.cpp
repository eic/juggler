// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

#include "EICDataSvc.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(EICDataSvc)

/// Standard Constructor
EICDataSvc::EICDataSvc(const std::string& name, ISvcLocator* svc) : PodioDataSvc(name, svc) {
  declareProperty("inputs", m_filenames = {}, "Names of the files to read");
  declareProperty("input", m_filename = "", "Name of the file to read");
}

/// Standard Destructor
EICDataSvc::~EICDataSvc() = default;
