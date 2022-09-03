#pragma once

#include <map>
#include <string>

#include <algorithms/property.h>

// Add boilerplate to service class definitions
//  - singleton --> no public constructor, no copy constructor, no assigmnent operator
//  - constructor should be protected so we can actually inherit from this class if needed
//    (mostly needed for the service base class)
#define ALGORITHMS_DEFINE_SERVICE(className) \
  protected: \
     className() : ServiceMixin<className>(#className) {} \
  public: \
     friend class ServiceMixin<className>; \
     className(const className&) = delete; \
     void operator=(const className&) = delete;
     
namespace algorithms {


// Service service --> keeps track of all services :)
//                 --> exposes the underlying property mixin of the service so
//                     we can configure the services by name in the framework 
//                     boundary plugin
//                 --> Special service (does not use the service mixin to avoid
//                     circularity)
class ServiceSvc {
  public:
    using ServiceMapType = std::map<std::string_view, PropertyMixin*>;
    static ServiceSvc& instance() {
      // This is guaranteed to be thread-safe from C++11 onwards.
      static ServiceSvc svc;
      return svc;
    }
    void add(std::string_view name, PropertyMixin* svcConfig) {
      m_services[name] = svcConfig;
    }
    PropertyMixin* get(std::string_view name) const {
      return m_services.at(name);
    }
    const ServiceMapType& get() const {
      return m_services;
    }

  private: 
     ServiceSvc() = default;
  public: 
     ServiceSvc(const ServiceSvc&) = delete;
     void operator=(const ServiceSvc&) = delete;
  private:
    ServiceMapType m_services;
};

// Thread-safe lazy-evaluated minimal service system
// CRTP mixin to add the instance method
// This could have been part of DEFINE_SERVICE macro, but I think it is better
// to keep the macro magic to a minimum to maximize transparency
template <class SvcType> class ServiceMixin : public PropertyMixin {
  public:
    static SvcType& instance() {
      // This is guaranteed to be thread-safe from C++11 onwards.
      static SvcType svc;
      return svc;
    }
    // constructor for the service mixin registers the service, except
    // for the ServiceSvc which is its own thing (and that would be ciricular)
    ServiceMixin(std::string_view name) {
      ServiceSvc::instance().add(name, this);
    }
};


}
