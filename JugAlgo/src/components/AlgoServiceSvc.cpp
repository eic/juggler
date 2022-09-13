// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten

#include <string>

#include <GaudiKernel/Service.h>
#include <JugAlgo/IAlgoServiceSvc.h>
#include <algorithms/logger.h>
#include <algorithms/service.h>

// too many services? :P
class AlgoServiceSvc : public extends<Service, IAlgoServiceSvc> {
public:
  AlgoServiceSvc(const std::string& name, ISvcLocator* svc) : base_class(name, svc) {}
  virtual ~AlgoServiceSvc() = default;

  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final { return StatusCode::SUCCESS; }
};

// Implementation

StatusCode AlgoServiceSvc::initialize() {
  StatusCode sc = Service::initialize();
  if (!sc.isSuccess()) {
    fatal() << "Error initializing ParticleSvc" << endmsg;
    return sc;
  }

  auto& serviceSvc = algorithms::ServiceSvc::instance();
  info() << "ServiceSvc declared " << serviceSvc.services().size() << " services" << endmsg;
  // loop over all services and handle each properly
  for (auto [name, svc] : serviceSvc.services()) {
    if (name == "LogSvc") {
      // TODO register action here
      info() << "LogSvc requested" << static_cast<unsigned>(static_cast<algorithms::LogSvc*>(svc)->defaultLevel())
             << endmsg;
    }
  }

  info() << "AlgoServiceSvc initialized successfully" << endmsg;
  return StatusCode::SUCCESS;
}
