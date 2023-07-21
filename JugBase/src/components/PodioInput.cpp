// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

#include "PodioInput.h"

#include "TFile.h"
#include "TROOT.h"

#include "JugBase/DataWrapper.h"
#include "JugBase/PodioDataSvc.h"

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(PodioInput)

PodioInput::PodioInput(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc), m_podioDataSvc(nullptr) {}

StatusCode PodioInput::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }

  // check whether we have the PodioDataSvc active
  m_podioDataSvc = dynamic_cast<PodioDataSvc*>(evtSvc().get());
  if (m_podioDataSvc == nullptr) {
    return StatusCode::FAILURE;
  }

  // TODO: add an upfront check for existence of data products

  return StatusCode::SUCCESS;
}

StatusCode PodioInput::execute() {
  [[maybe_unused]]
  size_t cntr = 0;
  // Re-create the collections from ROOT file

  for (auto& collName : m_collectionNames) {
    debug() << "Registering collection to read " << collName << endmsg;
    if (m_podioDataSvc->readCollection(collName).isFailure()) {
      return StatusCode::FAILURE;
    }
  }

  // Tell data service that we are done with requested collections
  m_podioDataSvc->endOfRead();
  return StatusCode::SUCCESS;
}

StatusCode PodioInput::finalize() {
  if (GaudiAlgorithm::finalize().isFailure()) {
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}
