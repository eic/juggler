#pragma once

#include <array>
#include <functional>
#include <ios>
#include <iostream>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string>

#include <algorithms/service.h>

// Simple thread-safe logger with optional overrides by the calling framework
namespace algorithms {

enum class LogLevel : unsigned { kJunk = 0, kDebug = 1, kInfo = 2, kWarning = 3, kError = 4 };
constexpr std::string_view logLevelName(LogLevel level) {
  // Compiler will warn if not all of the enum is covered
  switch (level) {
  case LogLevel::kJunk:
    return "JUNK";
  case LogLevel::kDebug:
    return "DEBUG";
  case LogLevel::kInfo:
    return "INFO";
  case LogLevel::kWarning:
    return "WARNING";
  case LogLevel::kError:
    return "ERROR";
  }
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
    // fmt::print("%s [%s] %s\n", logLevelName(l), caller, msg);
    std::cout << logLevelName(l) << " [" << caller << "] " << msg << std::endl;
  };
  ALGORITHMS_DEFINE_SERVICE(LogSvc)
};

namespace detail {
  // Output buffer that calls our global logger's report() function
  class LoggerBuffer : public std::stringbuf {
  public:
    LoggerBuffer(const LogLevel l, std::string_view caller)
        : m_mylevel{l}, m_caller{caller}, m_logger{algorithms::LogSvc::instance()} {}
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
      : m_caller{caller} {
    level(threshold);
  }

public:
  // Not done through Properties, as that is the responsible of the base Algo or Service
  void level(const LogLevel threshold) {
    m_level = threshold;
    m_error.threshold(m_level);
    m_warning.threshold(m_level);
    m_info.threshold(m_level);
    m_debug.threshold(m_level);
    m_junk.threshold(m_level);
  }
  LogLevel level() const { return m_level; }

protected:
  detail::LoggerStream& error() const { return m_error; }
  detail::LoggerStream& warning() const { return m_warning; }
  detail::LoggerStream& info() const { return m_info; }
  detail::LoggerStream& debug() const { return m_debug; }
  detail::LoggerStream& junk() const { return m_junk; }

private:
  const std::string m_caller;
  LogLevel m_level;
  mutable detail::LoggerStream m_error{m_caller, LogLevel::kError};
  mutable detail::LoggerStream m_warning{m_caller, LogLevel::kWarning};
  mutable detail::LoggerStream m_info{m_caller, LogLevel::kInfo};
  mutable detail::LoggerStream m_debug{m_caller, LogLevel::kDebug};
  mutable detail::LoggerStream m_junk{m_caller, LogLevel::kJunk};
};

} // namespace algorithms

#define endmsg std::flush
