#include "TestIRT.h"
#include "GaudiKernel/MsgStream.h"

namespace Jug::PID {
  TestIRT::TestIRT(const std::string& name, ISvcLocator* ploc)
      : GaudiAlgorithm(name, ploc) {}

  StatusCode TestIRT::initialize() {

    StatusCode sc = Algorithm::initialize();
    if (sc.isFailure())
      return sc;

    info() << "TestIRT: Inilializing..." << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode TestIRT::execute() {
    info() << "TestIRT: executing..." << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode TestIRT::finalize() {
    info() << "TestIRT: Finalizing..." << endmsg;
    return Algorithm::finalize(); // must be executed last
  }

  DECLARE_COMPONENT(TestIRT)
} // namespace Jug::PID
