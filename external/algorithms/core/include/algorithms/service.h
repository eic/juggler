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
  void add(std::string_view name, ServiceBase* svc) { m_services[name] = svc; }
  template <class T = ServiceBase> T* service(std::string_view name) const {
    return static_cast<T*>(m_services.at(name));
  }
  const ServiceMapType& services() const { return m_services; }
  bool active(std::string_view name) const { return m_services.count(name); }
  // Check if all service are in the `ready` state and all properties are set
  // TODO FIXME move to implementation file
  void validate() const {
    std::vector<std::string_view> not_ready;
    std::map<std::string_view, std::vector<std::string_view>> missing_props;
    for (const auto& [name, svc] : m_services) {
      if (!svc->ready()) {
        not_ready.push_back(name);
      }
      auto missing = svc->missingProperties();
      if (!missing.empty()) {
        missing_props[name] = missing;
      }
    }
    std::string err = "";
    if (!not_ready.empty()) {
      err += fmt::format("Encountered unitialized services: {}\n", not_ready);
    }
    if (!missing_props.empty()) {
      err += fmt::format("Encountered missing service properties: {}\n", missing_props);
    }
    if (!err.empty()) {
      throw ServiceError(fmt::format("Error initializing all services:\n{}", err));
    }
  }

private:
  ServiceSvc() : NameMixin{"ServiceSvc", "Special service that keeps track of all services"} {}

public:
  ServiceSvc(const ServiceSvc&) = delete;
  void operator=(const ServiceSvc&) = delete;

private:
  ServiceMapType m_services;
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
  Service(std::string_view name) : NameMixin{name, name} { ServiceSvc::instance().add(name, this); }
};

} // namespace algorithms

