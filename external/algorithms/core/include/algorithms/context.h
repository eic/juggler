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
  Context& operator=(const Context& rhs) = default;

  constexpr Context(const run_type r, const event_type e, const time_type t, bool is_mc = false)
      : m_run{r}, m_event{e}, m_time{t}, m_mc{is_mc} {}
  constexpr Context(const run_type r, const event_type e, bool is_mc = false)
      : m_run{r}, m_event{e}, m_mc{is_mc} {}
  constexpr Context(const time_type t) : m_time{t} {}

  void set(run_type r, event_type e, time_type t, bool is_mc) { *this = Context(r, e, t, is_mc); }
  void run(const run_type r) { m_run = r; }
  void event(const event_type e) { m_event = e; }
  void time(const time_type t) { m_time = t; }
  void mc(const bool b) { m_mc = b; }

  constexpr auto run() const { return m_run; }
  constexpr auto event() const { return m_event; }
  constexpr auto time() const { return m_time; }
  constexpr bool mc() const { return m_mc; }

  constexpr bool validEvent() const { return (m_run != kInvalidRun && m_event != kInvalidEvent); }
  constexpr bool validTime() const { return m_time != kInvalidTime; }
  constexpr bool valid() const { return validEvent() || validTime(); }

  // Get a context-specific ID
  constexpr uint64_t id() const {
    if (validEvent()) {
      return static_cast<uint64_t>(m_run) << 32 | m_event;
    } else if (validTime()) {
      return m_time;
    } else {
      throw ContextError("Cannot get uniqueID for invalid context. A valid context has at least "
                         "either a valid run+event, or a valid time stamp.");
    }
  }

  // framework-specific extension
  template <class T> void extension(T&& ext) { m_extension = std::forward<T>(ext); }
  template <class T> const auto& extension() const {
    return std::any_cast<const std::decay_t<T>&>(m_extension);
  }
  bool hasExtension() const { return m_extension.has_value(); }

private:
  run_type m_run     = kInvalidRun;
  event_type m_event = kInvalidEvent;
  time_type m_time   = kInvalidTime;
  bool m_mc          = false;

  std::any m_extension; // Framework dependent context
};

} // namespace algorithms
