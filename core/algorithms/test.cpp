#include <algorithms/logger.h>

using namespace algorithms;

// TODO implement PropertyBase? or use boost::any? need something better for setProperty
// so we don't need to repeat the type unnecessarily. Hard without type erasure though...

class testLogger : public LoggerMixin {
  public:
  testLogger() 
    : LoggerMixin("testLogger")
  {
    std::cout << "should display an info message" << std::endl;
    info() << "1" << endmsg;
    std::cout << "next debug message should not display" << std::endl;
    debug() << "2" << endmsg;
    level(LogLevel::kJunk);
    std::cout << "Changed the log level to junk, so the following message should display" << std::endl;
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
      info() << "integer property: " << getProperty<int>("integer") << endmsg;
      info() << "foo property: " << getProperty<std::string>("foo") << endmsg;
      if (hasProperty("bar")) {
        info() << "bar property: " << getProperty<std::string>("bar") << endmsg;
      }
      if (hasProperty("double")) {
        info() << "double (non-def.) property: " << getProperty<double>("double") << endmsg;
      }
    }
  private:
    Property<int> m_int {this, "integer", 3};
    Property<std::string> m_foo{this, "foo", "foo_property_value"};
    // and one without a default
    Property<std::string> m_bar{this, "bar"};
    Property<double> m_double{this, "double"};

};

int main() {
  testLogger tl;

  std::cout << "test of property, will print set properties. Should only print 2 to start (no bar, as it does not have a default value\n";
  propTest tp;
  tp.print();
  tp.setProperty<int>("integer", 10);
  tp.print();
  try {
    tp.setProperty<int>("foo", 1.);
    std::cout << "Should not have been reached, as foo does not exist" << std::endl;
    return 1;
  } catch (...) {
    std::cout << "setting a property that didn't exist failed popperly" << std::endl;
  } 
  std::cout << "Let's now also set a value for the third and fourth property so it gets printed as well" << std::endl;
  tp.setProperty<std::string>("bar", "bar_value");
  tp.setProperty<double>("double", 3.1415f);
  tp.print();


  return 0;
}
