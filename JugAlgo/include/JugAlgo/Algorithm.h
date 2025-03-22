// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten

#include <string>
#include <type_traits>

#include <algorithms/algorithm.h>
#include <algorithms/type_traits.h>

#include <Gaudi/Algorithm.h>
#include <GaudiKernel/Service.h>
#include <JugAlgo/IAlgoServiceSvc.h>
#include <JugAlgo/detail/DataProxy.h>

#define JUGALGO_DEFINE_ALGORITHM(algorithm, impl, ns)                                              \
  using algorithm##Base = Jug::Algo::Algorithm<impl>;                                              \
  namespace ns {                                                                                   \
  struct algorithm : algorithm##Base {                                                             \
    algorithm(const std::string& name, ISvcLocator* svcLoc) : algorithm##Base(name, svcLoc) {}     \
  };                                                                                               \
  DECLARE_COMPONENT(algorithm)                                                                     \
  }

namespace Jug::Algo {

template <class AlgoImpl> class Algorithm : public Gaudi::Algorithm {
public:
  using algo_type   = AlgoImpl;
  using input_type  = typename algo_type::input_type;
  using output_type = typename algo_type::output_type;
  using Input       = typename algo_type::Input;
  using Output      = typename algo_type::Output;

  Algorithm(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc)
      , m_algo{name}
      , m_output{this, m_algo.outputNames()}
      , m_input{this, m_algo.inputNames()} {
    defineProperties();
  }

  StatusCode initialize() override {
    debug() << "Initializing " << name() << endmsg;

    // Algorithms uses exceptions, Gaudi uses StatusCode --> catch and propagate
    try {
      // Grab the AlgoServiceSvc
      m_algo_svc = service("AlgoServiceSvc");
      if (!m_algo_svc) {
        error() << "Unable to get an instance of the AlgoServiceSvc" << endmsg;
        return StatusCode::FAILURE;
      }

      // Forward the log level of this algorithm
      const algorithms::LogLevel level{
          static_cast<algorithms::LogLevel>(msgLevel() > 0 ? msgLevel() - 1 : 0)};
      debug() << "Setting the logger level to " << algorithms::logLevelName(level) << endmsg;
      m_algo.level(level);

      // Init our data structures
      debug() << "Initializing data structures" << endmsg;
      m_input.init();
      m_output.init();

      // configure properties
      debug() << "Configuring properties" << endmsg;
      initProperties();

      // validate properties
      debug() << "Validating properties" << endmsg;
      m_algo.validate();

      // call the internal algorithm init
      debug() << "Initializing underlying algorithm " << m_algo.name() << endmsg;
      m_algo.init();
    } catch (const std::exception& e) {
      fatal() << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    try {
      m_algo.process(m_input.get(), m_output.get());
    } catch (const std::exception& e) {
      error() << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

protected:
  template <typename T> void setAlgoProp(std::string_view name, T&& value) {
    m_algo.template setProperty<T>(name, value);
  }
  template <typename T> T getAlgoProp(std::string name) const {
    return m_algo.template getProperty<T>(name);
  }
  bool hasAlgoProp(std::string_view name) const { return m_algo.hasProperty(name); }

private:
  // to be called during construction (before init())
  void defineProperties() {
    for (const auto& [key, prop] : m_algo.getProperties()) {
      std::visit(
          [this, key = key](auto&& val) {
            using T = std::decay_t<decltype(val)>;
            this->m_props.emplace(
                key, std::make_unique<Gaudi::Property<T>>(this, std::string(key), val));
          },
          prop.get());
    }
  }

  // to be called during init() --> will actually set the underlying alo properties
  void initProperties() {
    for (const auto& [key, prop] : m_algo.getProperties()) {
      std::visit(
          [this, key = key](auto&& val) {
            using T                = std::decay_t<decltype(val)>;
            const auto* gaudi_prop = static_cast<Gaudi::Property<T>*>(this->m_props.at(key).get());
            const auto prop_val    = gaudi_prop->value();
            this->m_algo.setProperty(key, prop_val);
          },
          prop.get());
    }
  }

  algo_type m_algo;
  SmartIF<IAlgoServiceSvc> m_algo_svc;
  detail::DataProxy<output_type> m_output;
  detail::DataProxy<input_type> m_input;
  std::map<std::string_view, std::unique_ptr<Gaudi::Details::PropertyBase>> m_props;
};

} // namespace Jug::Algo

