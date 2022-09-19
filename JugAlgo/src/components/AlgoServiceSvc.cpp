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

DECLARE_COMPONENT(AlgoServiceSvc)

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
  // Note: this code is kind of dangerous, as getting the types wrong will lead to
  // undefined runtime behavior.
  for (auto [name, svc] : serviceSvc.services()) {
    if (name == "LogSvc") {
      auto* logger = static_cast<algorithms::LogSvc*>(svc);
      const algorithms::LogLevel level{static_cast<algorithms::LogLevel>(msgLevel() > 0 ? msgLevel() - 1 : 0)};
      info() << "Setting up algorithms::LogSvc with default level " << algorithms::logLevelName(level) << endmsg;
      logger->defaultLevel(level);
      logger->action([this](const algorithms::LogLevel l, std::string_view caller, std::string_view msg) {
        const std::string text = fmt::format("[{}] {}", caller, msg);
        if (l == algorithms::LogLevel::kCritical) {
          this->fatal() << text << endmsg;
        } else if (l == algorithms::LogLevel::kError) {
          this->error() << text << endmsg;
        } else if (l == algorithms::LogLevel::kWarning) {
          this->warning() << text << endmsg;
        } else if (l == algorithms::LogLevel::kInfo) {
          this->info() << text << endmsg;
        } else if (l == algorithms::LogLevel::kDebug) {
          this->debug() << text << endmsg;
        } else if (l == algorithms::LogLevel::kTrace) {
          this->verbose() << text << endmsg;
        }
      });
      // set own log level to verbose so we actually display everything that is requested
      // (this overrides what was initally set through the OutputLevel property)
      updateMsgStreamOutputLevel(MSG::VERBOSE);
    }
  }

  info() << "AlgoServiceSvc initialized successfully" << endmsg;
  return StatusCode::SUCCESS;
}
