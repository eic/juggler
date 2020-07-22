#pragma once

#include "GaudiAlg/GaudiAlgorithm.h"

class HelloWorld : public GaudiAlgorithm {
public:
  HelloWorld(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize() override;
  StatusCode execute() override;
  StatusCode finalize() override;
};
