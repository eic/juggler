// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Logging service, which will use a callback to use the framework logging infrastructure.
// use ::action(void(LogLevel, std::string_view, std::string_view)) to register
// a logger.
//
// Also provides the LoggerMixin and LoggedService base classes
//
#pragma once

#include <array>
#include <functional>
#include <ios>
#include <iostream>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string>

#include <algorithms/error.h>
#include <algorithms/service.h>
#include <fmt/format.h>

// Simple thread-safe logger with optional overrides by the calling framework
namespace algorithms {

enum class LogLevel : unsigned {
  kTrace    = 0,
  kDebug    = 1,
  kInfo     = 2,
  kWarning  = 3,
  kError    = 4,
  kCritical = 5,
};
constexpr std::string_view logLevelName(LogLevel level) {
  // Compiler can warn if not all of the enum is covered
  switch (level) {
  case LogLevel::kTrace:
    return "TRACE";
  case LogLevel::kDebug:
    return "DEBUG";
  case LogLevel::kInfo:
    return "INFO";
  case LogLevel::kWarning:
    return "WARNING";
  case LogLevel::kError:
    return "ERROR";
  case LogLevel::kCritical:
    return "CRITICAL";
  }
  // Default return to make gcc happy, will never happen
  return "UNKNOWN";
}

// Note: the log action is responsible for dealing with concurrent calls
//       the default LogAction is a thread-safe example
class LogSvc : public Service<LogSvc> {
public:
  using LogAction = std::function<void(LogLevel, std::string_view, std::string_view)>;
  void defaultLevel(const LogLevel l) { m_level.set(l); }
  LogLevel defaultLevel() const { return m_level; }
  void action(LogAction a) { m_action = a; }
  void report(const LogLevel l, std::string_view caller, std::string_view msg) const {
    m_action(l, caller, msg);
  }

private:
  Property<LogLevel> m_level{this, "defaultLevel", LogLevel::kInfo};
  LogAction m_action = [](const LogLevel l, std::string_view caller, std::string_view msg) {
    static std::mutex m;
    std::lock_guard<std::mutex> lock(m);
    fmt::print("{} [{}] {}\n", logLevelName(l), caller, msg);
  };
  ALGORITHMS_DEFINE_SERVICE(LogSvc)
};

namespace detail {
  // Output buffer that calls our global logger's report() function
  class LoggerBuffer : public std::stringbuf {
  public:
    LoggerBuffer(const LogLevel l, std::string_view caller)
        : m_mylevel{l}, m_caller{caller}, m_logger{LogSvc::instance()} {}
    virtual int sync() {
      // report should deal with concurrency (the minimal version does)
      m_logger.report(m_mylevel, m_caller, this->str());
      this->str("");
      return 0;
    }

  private:
    // The output buffer knows the log level of its associated logger
    // (eg. is this the debug logger?)
    LogLevel m_mylevel;
    const std::string m_caller;
    const LogSvc& m_logger;
  };
  // thread-safe output stream for the logger
  class LoggerStream {
  public:
    LoggerStream(std::string_view caller, const LogLevel level,
                 const LogLevel threshold = LogSvc::instance().defaultLevel())
        : m_buffer{level, caller}, m_os{&m_buffer}, m_level{level}, m_threshold{threshold} {}
    LoggerStream()                    = delete;
    LoggerStream(const LoggerStream&) = delete;

    template <class Arg> LoggerStream& operator<<(Arg&& streamable) {
      if (m_level >= m_threshold) {
        std::lock_guard<std::mutex> lock{m_mutex};
        m_os << std::forward<Arg>(streamable);
        return *this;
      }
      return *this;
    }
    // To support input manipulators such as std::endl
    // Note: would be better with Concepts
    using IOManipType1 = std::ostream&(std::ostream&); // this capturs std::endl;
    using IOManipType2 = std::ios_base&(std::ios_base&);
    LoggerStream& operator<<(IOManipType1* f) {
      if (m_level >= m_threshold) {
        std::lock_guard<std::mutex> lock{m_mutex};
        f(m_os);
        return *this;
      }
      return *this;
    }
    LoggerStream& operator<<(IOManipType2* f) {
      if (m_level >= m_threshold) {
        std::lock_guard<std::mutex> lock{m_mutex};
        f(m_os);
        return *this;
      }
      return *this;
    }
    LogLevel threshold() const { return m_threshold; }
    void threshold(const LogLevel th) { m_threshold = th; }

  private:
    std::mutex m_mutex;
    LoggerBuffer m_buffer;
    std::ostream m_os;
    const LogLevel m_level;
    LogLevel m_threshold;
  };
} // namespace detail

// Mixin meant to add utility logger functions to algorithms/services/etc
class LoggerMixin {
public:
  LoggerMixin(std::string_view caller, const LogLevel threshold = LogSvc::instance().defaultLevel())
      : m_caller{caller}, m_logger{LogSvc::instance()} {
    level(threshold);
  }

public:
  // Not done through Properties, as that would require entanglement with the
  // PropertyMixin which is not appropriate here. It's the FW responsibility to set this
  // on the algorithm level if desired, before or during the init() stage.
  void level(const LogLevel threshold) {
    m_level = threshold;
    m_critical.threshold(m_level);
    m_error.threshold(m_level);
    m_warning.threshold(m_level);
    m_info.threshold(m_level);
    m_debug.threshold(m_level);
    m_trace.threshold(m_level);
  }
  LogLevel level() const { return m_level; }

protected:
  detail::LoggerStream& critical() const { return m_critical; }
  detail::LoggerStream& error() const { return m_error; }
  detail::LoggerStream& warning() const { return m_warning; }
  detail::LoggerStream& info() const { return m_info; }
  detail::LoggerStream& debug() const { return m_debug; }
  detail::LoggerStream& trace() const { return m_trace; }

  void critical(std::string_view msg) { report<LogLevel::kCritical>(msg); }
  void error(std::string_view msg) { report<LogLevel::kError>(msg); }
  void warning(std::string_view msg) { report<LogLevel::kWarning>(msg); }
  void info(std::string_view msg) { report<LogLevel::kInfo>(msg); }
  void debug(std::string_view msg) { report<LogLevel::kDebug>(msg); }
  void trace(std::string_view msg) { report<LogLevel::kTrace>(msg); }

  bool aboveCriticalThreshold() const { return m_level >= LogLevel::kCritical; }
  bool aboveErrorThreshold() const { return m_level >= LogLevel::kError; }
  bool aboveWarningThreshold() const { return m_level >= LogLevel::kWarning; }
  bool aboveInfoThreshold() const { return m_level >= LogLevel::kInfo; }
  bool aboveDebugThreshold() const { return m_level >= LogLevel::kDebug; }
  bool aboveTraceThreshold() const { return m_level >= LogLevel::kTrace; }

  // LoggerMixin also provides nice error raising
  // ErrorTypes needs to derive from Error, and needs to have a constructor that takes two
  // strings --> TODO add C++20 Concept version
  template <class ErrorType = Error> void raise(std::string_view msg) const {
    error() << msg << endmsg;
    throw ErrorType(msg);
  }

  // Endmsg (flush) support to initiate a sync() action
  static std::ostream& endmsg(std::ostream& os) {
    os.flush();
    return os;
  }

private:
  template <LogLevel l> void report(std::string_view msg) {
    if (l >= m_level) {
      m_logger.report(l, m_caller, msg);
    }
  }

  const std::string m_caller;
  LogLevel m_level;
  mutable detail::LoggerStream m_critical{m_caller, LogLevel::kCritical};
  mutable detail::LoggerStream m_error{m_caller, LogLevel::kError};
  mutable detail::LoggerStream m_warning{m_caller, LogLevel::kWarning};
  mutable detail::LoggerStream m_info{m_caller, LogLevel::kInfo};
  mutable detail::LoggerStream m_debug{m_caller, LogLevel::kDebug};
  mutable detail::LoggerStream m_trace{m_caller, LogLevel::kTrace};

  const LogSvc& m_logger;
};

// A Service base class with logger functionality
template <class SvcType> class LoggedService : public Service<SvcType>, public LoggerMixin {
protected:
  LoggedService(std::string_view name) : Service<SvcType>(name), LoggerMixin(name) {}
};

} // namespace algorithms

#define ALGORITHMS_DEFINE_LOGGED_SERVICE(className)                                                \
protected:                                                                                         \
  className() : LoggedService<className>(#className) {}                                            \
                                                                                                   \
public:                                                                                            \
  friend class Service<className>;                                                                 \
  className(const className&) = delete;                                                            \
  void operator=(const className&)   = delete;                                                     \
  constexpr static const char* kName = #className;

//#define endmsg std::flush
