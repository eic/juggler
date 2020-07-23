#include "HelloWorld.hh"
#include "GaudiKernel/MsgStream.h"

DECLARE_COMPONENT(HelloWorld)

HelloWorld::HelloWorld(const std::string& name, ISvcLocator* ploc)
    : GaudiAlgorithm(name, ploc) {}

StatusCode HelloWorld::initialize() {

  StatusCode sc = Algorithm::initialize();
  if (sc.isFailure())
    return sc;

  info() << "Hello World: Inilializing..." << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode HelloWorld::execute() {
  info() << "Hello World: ecuting..." << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode HelloWorld::finalize() {
  info() << "Hello World: Finalizing..." << endmsg;
  return Algorithm::finalize(); // must be executed last
}
