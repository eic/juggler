// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong

#include "JugBase/ACTSLogger.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/MsgStream.h"

namespace Acts {
std::unique_ptr<const Logger> getDefaultLogger(const std::string& name, const Logging::Level& lvl, std::ostream* /* unused */) {
  using namespace Logging;
  ServiceHandle<IMessageSvc> msgSvc("MessageSvc", name);
  msgSvc->setOutputLevel(lvl + 1);
  auto printPol = std::make_unique<GaudiPrintPolicy>(&(*msgSvc));
  printPol->setName(name);
  return std::make_unique<Acts::Logger>(std::move(printPol),
                                        std::make_unique<GaudiFilterPolicy>(&(*msgSvc), lvl));
}
}
