// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Sylvester Joosten
//
// Dummy that instantiates some objects to trigger compile errors early on in case of
// bugs. This should be migrated to tests TODO.
#include <algorithms/logger.h>

#include <algorithms/algorithm.h>
#include <algorithms/random.h>

#include <fmt/ranges.h>

using namespace algorithms;

class testLogger : public LoggerMixin {
public:
  testLogger() : LoggerMixin("testLogger") {
    std::cout << "should display an info message" << std::endl;
    info() << "1" << endmsg;
    std::cout << "next debug message should not display" << std::endl;
    debug() << "2" << endmsg;
    level(LogLevel::kTrace);
    std::cout << "Changed the log level to trace, so the following message should display"
              << std::endl;
    debug() << "3" << endmsg;
    std::cout << "error message should display" << std::endl;
    error() << "4" << endmsg;

    std::cout << std::endl;
    info() << "Checking the existing services:" << endmsg;
    for (const auto [key, value] : ServiceSvc::instance().services()) {
      info() << " - " << key << endmsg;
    }
  }
};

class propTest : public PropertyMixin, LoggerMixin {
public:
  propTest() : LoggerMixin("propTest") {}
  void print() const {
    info() << "integer property: " << m_int << endmsg;
    info() << "foo property: " << m_foo << endmsg;
    if (hasProperty("bar")) {
      info() << "bar property: " << m_bar << endmsg;
    }
    if (hasProperty("double")) {
      info() << "double (non-def.) property: " << m_double << endmsg;
    }
  }

private:
  Property<int> m_int{this, "integer", 3, "test integer property"};
  Property<std::string> m_foo{this, "foo", "foo_property_value", "test foo property"};
  // and one without a default
  Property<std::string> m_bar{this, "bar", "test bar property without default"};
  Property<double> m_double{this, "double", "test double property without default"};
};

int dummy() {
  testLogger tl;

  std::cout << "test of property, will print set properties. Should only print 2 to start (no bar, "
               "as it does not have a default value\n";
  propTest tp;
  tp.print();
  tp.setProperty("integer", 10);
  tp.print();
  try {
    tp.setProperty("foo", 1.);
    std::cout << "Should not have been reached, as foo does not exist" << std::endl;
    return 1;
  } catch (...) {
    std::cout << "setting a property that didn't exist failed popperly" << std::endl;
  }
  std::cout
      << "Let's now also set a value for the third and fourth property so it gets printed as well"
      << std::endl;
  tp.setProperty("bar", "bar_value");
  tp.setProperty("double", 3.1415f);
  tp.print();

  Algorithm<Input<int, double>, Output<double, std::vector<double>>> a{
      "myAlgo", {"int", "double"}, {"moredouble", "variable"}, "Test algo"};
  fmt::print("Algo input: {}\n", a.inputNames());
  fmt::print("Algo output: {}\n", a.outputNames());
  // The following won't compile as the number of strings (variable names) and
  // input/output tuples don't match
  // Algorithm<Input<int, double>, Output<double, std::vector<double>>> a{
  //    "myAlgo", {"int", "double"}, {"moredouble", "variable", "bar"}};

  return 0;
}
