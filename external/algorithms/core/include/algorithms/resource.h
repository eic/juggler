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

  // Copy constructor for the ResourceMixin needs to update the addresses of the contained
  // resources to refer to the new copies.
  ResourceMixin(const ResourceMixin& rhs) : m_resources(rhs.m_resources.size(), nullptr) {
    for (size_t i = 0; i < m_resources.size(); ++i) {
      m_resources[i] = rhs.m_resources[i]->relocate(this);
  //    std::cout << m_resources[i] << " " << rhs.m_resources[i] << std::endl;
  //    m_resources[i]->context(&m_context);
    }
  }
  ResourceMixin& operator=(const ResourceMixin& rhs) = delete;

  void context(const Context& c, const NameMixin* s = nullptr) {
    m_context = c;
    if (s) {
      m_context.scope(s);
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
    ResourceHandle(const ResourceMixin* owner)
        : m_offset{reinterpret_cast<const char*>(this) - reinterpret_cast<const char*>(owner)} {
      std::cout << "resources handle offset: " << m_offset << std::endl;
    }
    virtual void context(const Context*)   = 0;
    virtual const Context* context() const = 0;
    // return the relocated address for the copied object in the new instance
    ResourceHandle* relocate(ResourceMixin* clone) const {
      return reinterpret_cast<ResourceHandle*>(reinterpret_cast<char*>(clone) + m_offset);
    }

  private:
    // offset versus the ResourceMixin `this` pointer to allow relocating views
    const ptrdiff_t m_offset;
  };
  template <class ResourceType> class Resource : public ResourceHandle {
  public:
    template <class... Args>
    Resource(ResourceMixin* owner, Args&&... args)
        : ResourceHandle{owner}, m_impl{std::forward<Args>(args)...} {
      std::cout << "Hi, I'm the original and I'm here: " << this << std::endl;
      if (owner) {
        owner->registerResource(this);
        m_impl.context(&(owner->context()));
      } else {
        throw ResourceError(
            fmt::format("Attempting to create Resource '{}' without valid owner", m_impl.name()));
      }
    }
    Resource(const Resource& rhs) : m_impl{rhs.m_impl} {
      std::cout << "Hi I'm the copy and I'm here: " << this << " compared to previous: " << &rhs
                << " offset: " << (&rhs - this) << std::endl;
    }
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
};

} // namespace algorithms

