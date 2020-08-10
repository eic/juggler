#ifndef TYPE_HPP
#define TYPE_HPP

#include <string>
#include <typeinfo>

//https://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname 

namespace jug::helpers {

std::string demangle(const char* name);

template <class T>
std::string type(const T& t) {

    return demangle(typeid(t).name());
}
}

#endif
