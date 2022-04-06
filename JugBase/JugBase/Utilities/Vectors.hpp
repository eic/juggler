// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#ifndef JUGBASE_VectorHelpers_HH
#define JUGBASE_VectorHelpers_HH


#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Objects.h"


namespace Jug::Helpers {

  // Need some time to think about these...
  //template<typename T>
  //struct VectorToActs : public T {
  //  VectorToActs(const T& v) : T(v) {}
  //  operator Acts::Vector3() const {return {this->x(), this->y(), this->z()};}
  //};

  //template<typename T>
  //struct ArrayToRoot : public  T {
  //  ArrayToRoot(const T& v) : T(v) {}
  //  operator ROOT::Math::XYZVector() const {return {(*this)[0], (*this)[1], (*this)[1]}; }
  //};

}


#endif
