// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck, Benedikt Hegner

#include "PodioOutput.h"
#include "podio/podioVersion.h"
#include "GaudiKernel/ISvcLocator.h"
#include "JugBase/PodioDataSvc.h"
#include "rootUtils.h"

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(PodioOutput)

PodioOutput::PodioOutput(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc), m_firstEvent(true) {}

StatusCode PodioOutput::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }

  // check whether we have the PodioEvtSvc active
  m_podioDataSvc = dynamic_cast<PodioDataSvc*>(evtSvc().get());
  if (m_podioDataSvc == nullptr) {
    error() << "Failed to get the DataSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  m_framewriter = std::make_unique<podio::ROOTFrameWriter>(m_filename);
  m_switch       = KeepDropSwitch(m_outputCommands);
  return StatusCode::SUCCESS;
}

StatusCode PodioOutput::execute() {
  auto& frame = m_podioDataSvc->getEventFrame();

  // register for writing
  if (m_firstEvent) {
    auto collections = frame.getAvailableCollections();
    for (auto& collection_name : collections) {
      if (m_switch.isOn(collection_name)) {
	m_collection_names_to_write.push_back(collection_name);
      }
    }
    m_framewriter->writeFrame(frame, "events", m_collection_names_to_write);
  } else {
    m_framewriter->writeFrame(frame, "events", m_collection_names_to_write);
  }
  m_firstEvent = false;

  return StatusCode::SUCCESS;
}

/** PodioOutput::finalize
 * has to happen after all algorithms that touch the data store finish.
 * Here the job options are retrieved and stored to disk as a branch
 * in the metadata tree.
 *
 */
StatusCode PodioOutput::finalize() {
  info() << "Finalizing output algorithm" << endmsg;
  if (GaudiAlgorithm::finalize().isFailure()) {
    return StatusCode::FAILURE;
  }
  //// prepare job options metadata ///////////////////////
  // retrieve the configuration of the job
  // and write it to file as vector of strings
  debug() << "Preparing job options metadata" << endmsg;
  std::vector<std::string> config_data;
  const auto& jobOptionsSvc         = Gaudi::svcLocator()->getOptsSvc();
  const auto& configured_properties = jobOptionsSvc.items();
  for (const auto& per_property : configured_properties) {
    std::stringstream config_stream;
    // sample output:
    // HepMCToEDMConverter.genparticles = "GenParticles";
    // Note that quotes are added to all property values,
    // which leads to problems with ints, lists, dicts and bools.
    // For theses types, the quotes must be removed in postprocessing.
    config_stream << std::get<0>(per_property) << " = \"" << std::get<1>(per_property) << "\";" << std::endl;
    config_data.push_back(config_stream.str());
  }
  // Some default components are not captured by the job option service
  // and have to be traversed like this. Note that Gaudi!577 will improve this.
  debug() << "Appending default component metadata" << endmsg;
  for (const auto* name : {"ApplicationMgr", "MessageSvc", "NTupleSvc"}) {
    std::stringstream config_stream;
    auto svc = service<IProperty>(name);
    if (!svc.isValid()) {
      continue;
    }
    for (const auto* property : svc->getProperties()) {
      config_stream << name << "." << property->name() << " = \"" << property->toString() << "\";" << std::endl;
    }
    config_data.push_back(config_stream.str());
  }

  // Collect all the metadata
  podio::Frame config_metadata_frame{};
  config_metadata_frame.putParameter("gaudiConfigOptions", config_data);

  m_framewriter->writeFrame(config_metadata_frame, "configuration_metadata");

  auto& metadata_frame = m_podioDataSvc->getMetaDataFrame();
  m_framewriter->writeFrame(metadata_frame, "metadata");

  // write information into file
  m_framewriter->finish();

  return StatusCode::SUCCESS;
}
