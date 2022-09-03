#pragma once

// Add boilerplate to service class definitions
//  - singleton --> no public constructor, no copy constructor, no assigmnent operator
//  - constructor should be protected so we can actually inherit from this class if needed
//    (mostly needed for the service base class)
#define ALGORITHMS_DEFINE_SERVICE(className) \
  protected: \
     className() {} \
  public: \
     friend class ServiceMixin<className>; \
     className(const className&) = delete; \
     void operator=(const className&) = delete;
     
namespace algorithms {

// Thread-safe lazy-evaluated minimal service system
// CRTP mixin to add the instance method
// This could have been part of DEFINE_SERVICE macro, but I think it is better
// to keep the macro magic to a minimum to maximize transparency
template <class SvcType> class ServiceMixin {
  public:
    static SvcType& instance() {
      // This is guaranteed to be thread-safe from C++11 onwards.
      static SvcType svc;
      return svc;
    }
};

}
