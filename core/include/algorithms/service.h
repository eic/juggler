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
  // add a service to the ServiceSvc. Two options:
  // - prior to init --> just add the service to the service stack and set the initializer
  // - after init    --> add service to call stack, set the initializer,
  //                     manually call init for the service, revalidate
  template <class Svc>
  void add(
      ServiceBase* svc, const std::function<void(Svc&)>& init = [](Svc& s) { s.init(); }) {
    m_keys.push_back(Svc::kName);
    m_services[Svc::kName] = svc;
    // only add initializer if not already present (e.g. if for some reason the framework did
    // the registration early)
    if (m_initializers.count(Svc::kName) == 0) {
      setInit(init);
    }
    // call init if we are already in an initialized state (this is a straggler)
    if (m_init) {
      initSingle<Svc>();
      validate();
    }
  }

  template <class Svc> bool has() const { return m_services.count(Svc::kName); }
  template <class Svc> void setInit(const std::function<void(Svc&)>& init) {
    m_initializers[Svc::kName] = [init]() { init(Svc::instance()); };
  }

  // Loop over all services in order and initialize one-by-one
  // Finalize by validating that all is well
  void init() {
    // ensure we only call init once
    if (m_init) {
      throw ServiceError("Cannot initialize services twice");
    }
    // Call init for all the services and mark as ready
    for (const auto& name : m_keys) {
      try {
        m_initializers.at(name)();
      } catch (const std::exception& e) {
        // we encountered an issue, stop here so validation fails
        break;
      }
      auto svc = m_services.at(name);
      // Ensure our init made sense -- cannot have missing properties at this stage
      if (svc->missingProperties().size() > 0) {
        break;
      }
      svc->ready(true);
    }

    // Validate all services in case we encountered issues so we get useful error
    // reporting
    validate();

    // Label initialization as complete
    m_init = true;
  }
  // Single service init (meant to be called for stragglers after general init)
  template <class Svc> void initSingle() {
    std::string_view name = Svc::kName;
    if (!m_init) {
      throw ServiceError(
          fmt::format("initSingle<{}>() should not be called before/instead of init()", name));
    }
    auto svc = m_services.at(name);
    // ensure we only call init once
    if (svc->ready()) {
      throw ServiceError(fmt::format("Cannot initialize service {} - already initialized", name));
    }
    // call init
    m_initializers.at(name)();
    // mark as ready
    svc->ready(true);
  }

  template <class Svc = ServiceBase> Svc* service(std::string_view name) const {
    return static_cast<Svc*>(m_services.at(name));
  }
  // Check if all service properties are set
  // TODO FIXME move to implementation file
  void validate() const {
    std::map<std::string_view, std::vector<std::string_view>> missing_props;
    std::vector<std::string_view> not_initialized;
    for (std::string_view name : m_keys) {
      const auto svc = m_services.at(name);
      auto missing   = svc->missingProperties();
      if (!missing.empty()) {
        missing_props[name] = missing;
      }
      if (!svc->ready()) {
        not_initialized.push_back(name);
      }
      std::string err = "";
      if (missing_props.size() > 0) {
        err += fmt::format("Encountered missing service properties: {}\n", missing_props);
      }
      if (not_initialized.size() > 0) {
        err += fmt::format("Encountered uninitialized services: {}\n", not_initialized);
      }
      if (err.size() > 0) {
        throw ServiceError(fmt::format("Error initializing all services:\n{}", err));
      }
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
  bool m_init = false; // did we initialize the services already?
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

