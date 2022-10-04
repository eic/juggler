// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten

#pragma once

#include <GaudiKernel/IService.h>

// Juggler bindings for all required algorithms services
// Will setup all necessary services (using the algorithms::ServiceSvc)

class GAUDI_API IAlgoServiceSvc : virtual public IService {
public:
  DeclareInterfaceID(IAlgoServiceSvc, 1, 0);
  virtual ~IAlgoServiceSvc() {}
  // No actual API needed, as this service will do all the magic behind the screens
};
