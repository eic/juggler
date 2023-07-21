// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Benedikt Hegner

#ifndef JUGBASE_EICLEGACYDATASVC_H
#define JUGBASE_EICLEGACYDATASVC_H

#include "JugBase/PodioLegacyDataSvc.h"

class EICLegacyDataSvc : public PodioLegacyDataSvc {
public:
  /// Standard Constructor
  EICLegacyDataSvc(const std::string& name, ISvcLocator* svc);

  /// Standard Destructor
  virtual ~EICLegacyDataSvc();
};
#endif  // JUGBASE_EICLEGACYDATASVC_H
