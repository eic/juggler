#pragma once

#include <any>
#include <limits>

#include <algorithms/error.h>

namespace algorithms {

class ContextError : public Error {
public:
  ContextError(std::string_view msg) : Error{msg, "algorithms::ContextError"} {}
};

class Context {
public:
  using run_type   = uint32_t;
  using event_type = uint32_t;
  using time_type  = uint64_t;

  constexpr static const run_type kInvalidRun     = std::numeric_limits<run_type>::max();
  constexpr static const event_type kInvalidEvent = std::numeric_limits<event_type>::max();
  constexpr static const time_type kInvalidTime   = std::numeric_limits<time_type>::max();

  constexpr Context()     = default;
  Context(const Context&) = default;
  // Operator= does _not_ copy over scope, as scope persists between events. Scope has to
  // be set manually
  Context& operator=(const Context& rhs) {
    m_run   = rhs.m_run;
    m_event = rhs.m_event;
    m_time  = rhs.m_time;
    m_mc    = rhs.m_mc;
    // Scope _not_ copied
    return *this;
  }

  constexpr Context(const run_type r, const event_type e, const time_type t, bool is_mc = false)
      : m_run{r}, m_event{e}, m_time{t}, m_mc{is_mc} {}
  constexpr Context(const run_type r, const event_type e, bool is_mc = false)
      : m_run{r}, m_event{e}, m_mc{is_mc} {}
  constexpr Context(const time_type t) : m_time{t} {}

  // set does _not_ set Scope, as Scope persists between event context changes
  void set(run_type r, event_type e, time_type t, bool is_mc) { *this = Context(r, e, t, is_mc); }
  void run(const run_type r) { m_run = r; }
  void event(const event_type e) { m_event = e; }
  void time(const time_type t) { m_time = t; }
  void mc(const bool b) { m_mc = b; }
  void scope(const NameMixin& s) { m_scope = &s; }

  constexpr auto run() const { return m_run; }
  constexpr auto event() const { return m_event; }
  constexpr auto time() const { return m_time; }
  constexpr bool mc() const { return m_mc; }
  constexpr const NameMixin* scope() const { return m_scope; }

  constexpr bool validEvent() const { return (m_run != kInvalidRun && m_event != kInvalidEvent); }
  constexpr bool validTime() const { return m_time != kInvalidTime; }
  constexpr bool validScope() const { return m_scope != nullptr; }
  constexpr bool valid() const { return validEvent() || validTime(); }

  // Get a context-specific ID
  constexpr uint64_t id() const {
    const auto scopeID = (m_scope ? m_scope->id() : 0);
    if (validEvent()) {
      return (static_cast<uint64_t>(m_run) << 32 | m_event) - scopeID;
    } else if (validTime()) {
      return m_time - scopeID;
    } else {
      throw ContextError("Cannot get uniqueID for invalid context. A valid context has at least "
                         "either a valid run+event, or a valid time stamp.");
    }
  }

  // Describe this context
  std::string describe() const {
    return fmt::format("run: {}, event: {}, time: {}, scope: {}{}",
                       m_run != kInvalidRun ? std::to_string(m_run) : "unknown",
                       m_event != kInvalidEvent ? std::to_string(m_event) : "unknown",
                       m_time != kInvalidTime ? std::to_string(m_time) : "unknown",
                       m_scope != nullptr ? m_scope->name() : "global", m_mc ? " (MC data)" : "");
  }

  // framework-specific extension
  template <class T> void extension(T&& ext) { m_extension = std::forward<T>(ext); }
  template <class T> const auto& extension() const {
    return std::any_cast<const std::decay_t<T>&>(m_extension);
  }
  bool hasExtension() const { return m_extension.has_value(); }

private:
  run_type m_run           = kInvalidRun;
  event_type m_event       = kInvalidEvent;
  time_type m_time         = kInvalidTime;
  const NameMixin* m_scope = nullptr;
  bool m_mc                = false;

  std::any m_extension; // Framework dependent context
};

} // namespace algorithms
