// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Basic Service infrastructure, implemented as lazy-evaluated singletons (thread-safe).
//
// Provides the special ServiceSvc (that provides access to all instantiated services as
// Configurable*), and the ServiceMixin (and related ALGORITHMS_DEFINE_SERVICE macro).
//
#pragma once

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <map>
#include <string>

#include <algorithms/error.h>
#include <algorithms/name.h>
#include <algorithms/property.h>

// Add boilerplate to service class definitions
//  - singleton --> no public constructor, no copy constructor, no assigmnent operator
//  - constructor should be protected so we can actually inherit from this class if needed
//    (mostly needed for the service base class)
#define ALGORITHMS_DEFINE_SERVICE(className)                                                       \
protected:                                                                                         \
  className() : Service<className>(#className) {}                                                  \
                                                                                                   \
public:                                                                                            \
  friend class Service<className>;                                                                 \
  className(const className&) = delete;                                                            \
  void operator=(const className&)   = delete;                                                     \
  constexpr static const char* kName = #className;

namespace algorithms {

class ServiceError : public Error {
public:
  ServiceError(std::string_view msg) : Error{msg, "algorithms::ServiceError"} {}
};

class ServiceBase : public PropertyMixin {
public:
  // Is this service ready? (active and initialized)
  bool ready() const { return m_ready; }
  void ready(bool state) { m_ready = state; }

private:
  bool m_ready{false};
};

// Service service --> keeps track of all services :)
//                 --> exposes the underlying ServiceBase* object of the service so
//                     we can configure the services by name in the framework
//                     boundary plugin
//                 --> Special service (does not use the Service base class to avoid
//                     circularity)
class ServiceSvc : public NameMixin {
public:
  using ServiceMapType = std::map<std::string_view, ServiceBase*>;
  static ServiceSvc& instance() {
    // This is guaranteed to be thread-safe from C++11 onwards.
    static ServiceSvc svc;
    return svc;
  }
  template <class Svc>
  void add(
      ServiceBase* svc, const std::function<void(Svc&)>& init = [](Svc& s) { s.init(); }) {
    m_keys.push_back(Svc::kName);
    m_services[Svc::kName] = svc;
    setInit(init);
  }

  template <class Svc> bool has() const { return m_services.count(Svc::kName); }
  template <class Svc> void setInit(const std::function<void(Svc&)>& init) {
    std::cerr << "DBGDBG - Setting init callback for " << Svc::kName << std::endl;
    m_initializers[Svc::kName] = [init]() { init(Svc::instance()); };
  }

  // Loop over all services in order and initialize one-by-one
  // Finalize by validating that all is well
  void init() {
    // Validate the services before we call init
    validate();
    // Call init for all the services and mark as ready
    std::cerr << fmt::format("DBGDBG - Calling init for all requested algorithms: {}", m_keys)
              << std::endl;
    for (const auto& name : m_keys) {
      std::cerr << "DBGDBG - Initializing: " << name << std::endl;
      m_initializers[name]();
      m_services[name]->ready(true);
    }
  }

  template <class Svc = ServiceBase> Svc* service(std::string_view name) const {
    return static_cast<Svc*>(m_services.at(name));
  }
  // Check if all service properties are set
  // TODO FIXME move to implementation file
  void validate() const {
    std::map<std::string_view, std::vector<std::string_view>> missing_props;
    for (const auto& [name, svc] : m_services) {
      auto missing = svc->missingProperties();
      if (!missing.empty()) {
        missing_props[name] = missing;
      }
    }
    std::string err = "";
    if (!missing_props.empty()) {
      err += fmt::format("Encountered missing service properties: {}\n", missing_props);
      throw ServiceError(fmt::format("Error initializing all services:\n{}", err));
    }
  }
  const auto& services() const { return m_services; }
  bool ready() const {
    for (const auto& key : m_keys) {
      if (!m_services.at(key)->ready()) {
        return false;
      }
    }
    return true;
  }

private:
  ServiceSvc() : NameMixin{"ServiceSvc", "Special service that keeps track of all services"} {}

public:
  ServiceSvc(const ServiceSvc&) = delete;
  void operator=(const ServiceSvc&) = delete;

private:
  std::vector<std::string_view> m_keys;                // Ordered list of service keys
  std::map<std::string_view, ServiceBase*> m_services; // Map of services for easier lookup
  std::map<std::string_view, std::function<void()>> m_initializers; // Init calls
};

// Thread-safe lazy-evaluated minimal service system
// CRTP base class to add the instance method
// This could have been part of DEFINE_SERVICE macro, but I think it is better
// to keep the macro magic to a minimum to maximize transparency
template <class SvcType> class Service : public ServiceBase, public NameMixin {
public:
  static SvcType& instance() {
    // This is guaranteed to be thread-safe from C++11 onwards.
    static SvcType svc;
    return svc;
  }
  // constructor for the service base class registers the service, except
  // for the ServiceSvc which is its own thing (avoid circularity)
  // (services don't currently have an attached description)
  Service(std::string_view name) : NameMixin{name, name} {
    ServiceSvc::instance().add<SvcType>(this);
  }
};

} // namespace algorithms

