#pragma once

#include <functional>
#include <string>

namespace algorithms {

  class PropertyBase;
  class ServiceBase;

  class StatusCode {
    public:
      enum Value : uint {
        SUCCESS,
        FAILURE
      };
      StatusCode() = default;
      constexpr StatusCode(Value value): m_value(value) { };
      constexpr operator Value() const { return m_value; }
      explicit operator bool() const = delete;        
      constexpr bool operator==(StatusCode a) const { return m_value == a.m_value; }
      constexpr bool operator!=(StatusCode a) const { return m_value != a.m_value; }
      constexpr bool isFailure() const { return m_value == FAILURE; }
    private:
      Value m_value;
  };

  class AlgorithmBase {
  private:
    std::vector<std::pair<std::string, PropertyBase*>> m_properties;
    std::vector<std::pair<std::string, ServiceBase*>> m_services;
  public:
    AlgorithmBase() = default;
    void registerProperty(PropertyBase* property, const std::string& name) {
      m_properties.push_back(std::make_pair(name, property));
    }
    void registerService(ServiceBase* service, const std::string& name) {
      m_services.push_back(std::make_pair(name, service));
    }

  };

  template<typename Out, typename In>
  class JugAlgorithm : public AlgorithmBase {
  private:
    std::string m_name;

    static std::function<std::ostream&()> m_debug;
    static std::function<std::ostream&()> m_info;
    static std::function<std::ostream&()> m_warning;
    static std::function<std::ostream&()> m_error;
    static std::function<std::ostream&()> m_endmsg;

    static void SetLoggers(
      std::function<std::ostream&()> debug,
      std::function<std::ostream&()> info,
      std::function<std::ostream&()> warning,
      std::function<std::ostream&()> error
    ) {
      m_debug = debug;
      m_info = info;
      m_warning = warning;
      m_error = error;
    }

  public:
    JugAlgorithm(const std::string& name = "")
    : m_name(name) {
    }

    virtual Out operator()(const In&) const = 0;

  protected:
    static std::ostream& debug() { return m_debug(); };
    static std::ostream& info() { return m_info(); };
    static std::ostream& warning() { return m_warning(); };
    static std::ostream& error() { return m_error(); };
    static std::ostream& endmsg() { return m_endmsg(); };

  };

} // namespace algorithms
