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
  void context(const Context& c) { m_context = c; }
  const Context& context() const { return m_context; }

private:
  const SvcType& m_service;
  Context m_context;
};

// Mixin for Resource handling at the algorithm (or service) level (mirroring the
// PropertyMixin)
class ResourceMixin {
public:
  class ResourceHandle;
  void resourceContext(const Context& c) {
    for (auto r : m_resources) {
      r->context(c);
    }
  }

private:
  void registerResource(ResourceHandle& resource) { m_resources.emplace_back(&resource); }
  std::vector<gsl::not_null<ResourceHandle*>> m_resources;

public:
  // Resource wrapper that acts as a handle to the Resource that automatically takes care of Context
  // management. Implementation is simular to Property
  class ResourceHandle {
  public:
    virtual void context(const Context& c) = 0;
  };
  template <class ResourceType> class Resource : public ResourceHandle {
  public:
    Resource(ResourceMixin* owner) {
      if (owner) {
        owner->registerResource(*this);
      } else {
        throw ResourceError(
            fmt::format("Attempting to create Resource '{}' without valid owner", m_impl.name()));
      }
    }
    void context(const Context& c) final { m_impl.context(c); }

    // Indirection operators to work with the underlying resource
    ResourceType* operator->() { return &m_impl; }
    const ResourceType* operator->() const { return &m_impl; }
    ResourceType& operator*() { return m_impl; }
    const ResourceType& operator*() const { return m_impl; }

  private:
    ResourceType m_impl;
  };
};

} // namespace algorithms

