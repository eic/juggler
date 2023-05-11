// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Context-aware Resources provided by Services
//
// Provides a ResourceMixin, similar to the PropertyMixin, to provide user-friendly
// resource integration within algorithms
//
#pragma once

#include <algorithms/context.h>
#include <algorithms/error.h>
#include <algorithms/name.h>

#include <gsl/gsl>
#include <vector>

namespace algorithms {

class ResourceError : public Error {
public:
  ResourceError(std::string_view msg) : Error{msg, "algorithms::ResourceError"} {}
};

// Service Resource CRPT that keeps a handle to the service (ensuring that the service will be
// initialized). A Resource is aware of its context and communicates appropriately with the relevant
// Service
template <class SvcType> class SvcResource : public NameMixin {
public:
  SvcResource(std::string_view name) : NameMixin{name, name}, m_service{SvcType::instance()} {}
  const SvcType& service() const { return m_service; }
  void context(const Context* c) { m_context = c; }
  const Context* context() const { return m_context; }
  // Called whenever something changed (e.g. context update)
  // to be overriden
  void update() {}

private:
  const SvcType& m_service;
  const Context* m_context;
};

// Mixin for Resource handling at the algorithm (or service) level (mirroring the
// PropertyMixin)
class ResourceMixin {
public:
  class ResourceHandle;

  ResourceMixin() = default;

  // Copy constructor for the ResourceMixin assumes new auto-registration (as
  // pointers need to be relocated to the new instance)
  ResourceMixin(const ResourceMixin&) : m_resources() {}

  ResourceMixin& operator=(const ResourceMixin& rhs) = delete;

  void context(const Context& ctx, const NameMixin* scope = nullptr) {
    m_context = ctx;
    if (scope) {
      m_context.scope(scope);
    }
  }
  const Context& context() const { return m_context; }

private:
  void registerResource(gsl::not_null<ResourceHandle*> resource) {
    m_resources.emplace_back(resource);
  }
  std::vector<ResourceHandle*> m_resources;
  Context m_context;

public:
  // Resource wrapper that acts as a handle to the Resource that automatically takes care of Context
  // management. Implementation is similar to Property
  class ResourceHandle {
  public:
    ResourceHandle()                       = default;
    virtual void context(const Context*)   = 0;
    virtual const Context* context() const = 0;
  };
  template <class ResourceType> class Resource : public ResourceHandle {
  public:
    template <class... Args>
    Resource(ResourceMixin* owner, Args&&... args) : m_impl{std::forward<Args>(args)...} {
      std::cout << "DEBUGDEBUG Hi, I'm a resource and I'm here: " << this << std::endl;
      if (owner) {
        owner->registerResource(this);
        m_impl.context(&(owner->context()));
      } else {
        throw ResourceError(
            fmt::format("Attempting to create Resource '{}' without valid owner", m_impl.name()));
      }
    }
    Resource(const Resource& rhs) = delete;
    void context(const Context* c) final { m_impl.context(c); }
    const Context* context() const final { return m_impl.context(); }

    // Indirection operators to work with the underlying resource
    ResourceType* operator->() { return &m_impl; }
    const ResourceType* operator->() const { return &m_impl; }
    ResourceType& operator*() { return m_impl; }
    const ResourceType& operator*() const { return m_impl; }

  private:
    ResourceType m_impl;
  };
}; // namespace algorithms

} // namespace algorithms

