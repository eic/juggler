// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong
#ifndef JUGBASE_FCCDATASVC_H
#define JUGBASE_FCCDATASVC_H

#include "JugBase/PodioDataSvc.h"

class EICDataSvc : public PodioDataSvc {

public:
  /// Standard Constructor
  EICDataSvc(const std::string& name, ISvcLocator* svc);

  /// Standard Destructor
  virtual ~EICDataSvc();
};
#endif  // JUGBASE_FCCDATASVC_H
