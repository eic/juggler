// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten

#include <string>

#include <GaudiKernel/Service.h>
#include <JugAlgo/IAlgoServiceSvc.h>
#include <JugBase/IGeoSvc.h>

#include <algorithms/geo.h>
#include <algorithms/logger.h>
#include <algorithms/random.h>
#include <algorithms/service.h>

// too many services? :P
class AlgoServiceSvc : public extends<Service, IAlgoServiceSvc> {
public:
  AlgoServiceSvc(const std::string& name, ISvcLocator* svc) : base_class(name, svc) {}
  virtual ~AlgoServiceSvc() = default;

  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final { return StatusCode::SUCCESS; }

private:
  SmartIF<IGeoSvc> m_geoSvc;
  Gaudi::Property<size_t> m_randomSeed{this, "randomSeed", 1};
};

DECLARE_COMPONENT(AlgoServiceSvc)

// Implementation

StatusCode AlgoServiceSvc::initialize() {
  StatusCode sc = Service::initialize();
  if (!sc.isSuccess()) {
    fatal() << "Error initializing AlgoServiceSvc" << endmsg;
    return sc;
  }

  try {
    auto& serviceSvc = algorithms::ServiceSvc::instance();
    info() << "ServiceSvc declared " << serviceSvc.services().size() << " services" << endmsg;
    // Now register custom initializers
    // Starting with the logger
    const algorithms::LogLevel level{
        static_cast<algorithms::LogLevel>(msgLevel() > 0 ? msgLevel() - 1 : 0)};
    info() << "Setting up algorithms::LogSvc with default level " << algorithms::logLevelName(level)
           << endmsg;
    serviceSvc.setInit<algorithms::LogSvc>([this,level](auto&& logger) {
      this->info() << "Initializing the algorithms::LogSvc using the Gaudi logger" << endmsg;
      logger.defaultLevel(level);
      logger.init(
          [this](const algorithms::LogLevel l, std::string_view caller, std::string_view msg) {
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
    });
    // geo Service
    info() << "Setting up algorithms::GeoSvc" << endmsg;
    m_geoSvc = service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc in the right order in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }
    serviceSvc.setInit<algorithms::GeoSvc>([this](auto&& g) {
      this->info() << "Initializing algorithms::RandomSvc with the Juggler GeoSvc" << endmsg;
      g.init(m_geoSvc->detector());
    });
    // setup random service
    info() << "Setting up algorithms::RandomSvc\n"
           << "  --> using internal STL 64-bit MT engine\n"
           << "  --> seed set to" << m_randomSeed << endmsg;
    serviceSvc.setInit<algorithms::RandomSvc>([this](auto&& r) {
      this->info() << "Initializing the algorithms::RandomSvc" << endmsg;
      r.setProperty("seed", m_randomSeed);
      r.init();
    });
    info() << "Initializing services" << endmsg;

    // first set own log level to verbose so we actually display everything that is requested
    // (this overrides what was initally set through the OutputLevel property)
    updateMsgStreamOutputLevel(MSG::VERBOSE);

    // Validate our service setup (this will throw if we encounter issues)
    serviceSvc.init();

  } catch (const std::exception& e) {
    fatal() << e.what() << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "AlgoServiceSvc initialized successfully" << endmsg;
  return StatusCode::SUCCESS;
}
