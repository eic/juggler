#pragma once

#include <iostream>
#include <ostream>
#include <ios>
#include <sstream>
#include <string>
#include <functional>
#include <mutex>
#include <array>

#include <algorithms/service.h>
#include <fmt/core.h>

// Simple thread-safe logger with optional overrides by the calling framework
namespace algorithms {

enum class LogLevel : unsigned {
  kJunk = 0,
  kDebug = 1,
  kInfo = 2,
  kWarning = 3,
  kError = 4
};
constexpr std::array<const char*, 5> kLogLevelNames{
    "JUNK", "DEBUG", "INFO", "WARNING", "ERROR"};

// TODO: integrate proper properties to configure the default level
// Note: the log action is responsible for dealing with concurrent calls
//       the default LogAction is a thread-safe example
class LogSvc : public ServiceMixin<LogSvc> {
  public:
    using LogAction = std::function<void(LogLevel, std::string_view, std::string_view)>;
    void defaultLevel(const LogLevel l) {m_level = l;}
    LogLevel defaultLevel() const {return m_level;}
    void action(LogAction a) {
      m_action = a;
    }
    void report(const LogLevel l, std::string_view caller, std::string_view msg) const {
      m_action(l, caller, msg);
    }

  private:
    LogLevel m_level = LogLevel::kInfo;
    LogAction m_action = [](const LogLevel l, std::string_view caller, std::string_view msg) {
      static std::mutex m;
      std::lock_guard<std::mutex> lock(m);
      fmt::print("%s [%s] %s\n", kLogLevelNames[l], caller, msg);
      //std::cout << kLogLevelNames[static_cast<unsigned>(l)] << " [" << caller << "] " << msg << std::endl;
    };
  ALGORITHMS_DEFINE_SERVICE(LogSvc)
};

namespace detail {
  // Output buffer that calls our global logger's report() function
  class LoggerBuffer : public std::stringbuf {
    public:
      LoggerBuffer(const LogLevel l, std::string_view caller) : m_mylevel{l}, m_caller{caller}, m_logger{algorithms::LogSvc::instance()} {}
      virtual int sync() {
        // report should deal with concurrency (the minimal version does)
        m_logger.report(m_mylevel, m_caller, this->str());
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
      LoggerStream(std::string_view caller, const LogLevel level, const LogLevel threshold = LogSvc::instance().defaultLevel())
        : m_buffer{level, caller}
        , m_os{&m_buffer}
        , m_level{level}
        , m_threshold{threshold} {}
      LoggerStream() = delete;
      LoggerStream(const LoggerStream&) = delete;
      
      template <class Arg>
      LoggerStream& operator<<(Arg&& streamable) {
        if (m_level >= m_threshold) {
          std::lock_guard<std::mutex> lock{m_mutex};
          m_os << std::forward<Arg>(streamable);
          return *this;
        }
        return *this;
      }
      // To support input manipulators such as std::endl
      // Note: would be better with Concepts
      using IOManipType1 = std::ostream&(std::ostream&);  // this capturs std::endl;
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
      LogLevel threshold() const {return m_threshold;}
      void threshold(const LogLevel th) {m_threshold = th;}

    private:
      std::mutex m_mutex;
      LoggerBuffer m_buffer;
      std::ostream m_os;
      const LogLevel m_level;
      LogLevel m_threshold;
  };
}

// Mixin meant to add utility logger functions to algorithms/services/etc
// TODO: add property to configure log level instead of using member function
class LoggerMixin {
  public:
    LoggerMixin(std::string_view caller, const LogLevel threshold = LogSvc::instance().defaultLevel()) 
      : m_caller{caller} {
      level(threshold);
    }
  public:
      void level(const LogLevel threshold) { 
        m_level = threshold;
        m_error.threshold(m_level);
        m_warning.threshold(m_level);
        m_info.threshold(m_level);
        m_debug.threshold(m_level);
        m_junk.threshold(m_level);
      }
      LogLevel level() const {return m_level;}

  protected:
    detail::LoggerStream& error() {
      return m_error;
    }
    detail::LoggerStream& warning() {
      return m_warning;
    }
    detail::LoggerStream& info() {
      return m_info;
    }
    detail::LoggerStream& debug() {
      return m_debug;
    }
    detail::LoggerStream& junk() {
      return m_junk;
    }

  private:
    const std::string m_caller;
    LogLevel m_level;
    detail::LoggerStream m_error{m_caller, LogLevel::kError};
    detail::LoggerStream m_warning{m_caller, LogLevel::kWarning};
    detail::LoggerStream m_info{m_caller, LogLevel::kInfo};
    detail::LoggerStream m_debug{m_caller, LogLevel::kDebug};
    detail::LoggerStream m_junk{m_caller, LogLevel::kJunk};
};

}

#define endmsg std::flush
