#pragma once
#include "GaudiAlg/GaudiAlgorithm.h"
#include "JugBase/IGeoSvc.h"

namespace Jug::PID {
  class TestIRT : public GaudiAlgorithm {
    public:
      TestIRT(const std::string& name, ISvcLocator* pSvcLocator);
      StatusCode initialize() override;
      StatusCode execute() override;
      StatusCode finalize() override;
  };
} // namespace Jug::PID
