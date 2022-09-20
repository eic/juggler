#pragma once

#include <map>
#include <string>

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

// Service service --> keeps track of all services :)
//                 --> exposes the underlying Configurable object of the service so
//                     we can configure the services by name in the framework
//                     boundary plugin
//                 --> Special service (does not use the Service base class to avoid
//                     circularity)
class ServiceSvc {
public:
  using ServiceMapType = std::map<std::string_view, Configurable*>;
  static ServiceSvc& instance() {
    // This is guaranteed to be thread-safe from C++11 onwards.
    static ServiceSvc svc;
    return svc;
  }
  void add(std::string_view name, Configurable* svc) { m_services[name] = svc; }
  Configurable* service(std::string_view name) const { return m_services.at(name); }
  const ServiceMapType& services() const { return m_services; }
  bool active(std::string_view name) const { return m_services.count(name); }

private:
  ServiceSvc() = default;

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
template <class SvcType> class Service : public PropertyMixin {
public:
  static SvcType& instance() {
    // This is guaranteed to be thread-safe from C++11 onwards.
    static SvcType svc;
    return svc;
  }
  // constructor for the service base class registers the service, except
  // for the ServiceSvc which is its own thing (avoid circularity)
  Service(std::string_view name) : m_name{name} { ServiceSvc::instance().add(name, this); }
  std::string_view name() const { return m_name; }

private:
  const std::string m_name;
};

} // namespace algorithms

