// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck

#include "PodioOutput.h"
#include "podio/podioVersion.h"
#include "GaudiKernel/ISvcLocator.h"
#include "JugBase/PodioDataSvc.h"
#include "TFile.h"
#include "rootutils.h"

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(PodioOutput)

PodioOutput::PodioOutput(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc), m_firstEvent(true),
      m_podioDataSvc(nullptr), m_datatree(nullptr), m_metadatatree(nullptr),
      m_runMDtree(nullptr), m_evtMDtree(nullptr), m_colMDtree(nullptr) {}

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

  m_file = std::unique_ptr<TFile>(TFile::Open(m_filename.value().c_str(), "RECREATE", "data file"));
  // Both trees are written to the ROOT file and owned by it
  // PodioDataSvc has ownership of EventDataTree
  m_datatree = m_podioDataSvc->eventDataTree();
  m_datatree->SetDirectory(m_file.get());
  m_metadatatree = new TTree("metadata", "Metadata tree");
  m_runMDtree = new TTree("run_metadata", "Run metadata tree");
  m_evtMDtree = new TTree("evt_metadata", "Event metadata tree");
  m_colMDtree = new TTree("col_metadata", "Collection metadata tree");

  m_evtMDtree->Branch("evtMD", "GenericParameters", m_podioDataSvc->getProvider().eventMetaDataPtr() ) ;
  m_switch       = KeepDropSwitch(m_outputCommands);
  return StatusCode::SUCCESS;
}

void PodioOutput::resetBranches(const std::vector<std::pair<std::string, podio::CollectionBase*>>& collections) {
  for (const auto& [collName, collBuffers] : collections) {
    auto buffers = collBuffers->getBuffers();
    auto* data = buffers.data;
    auto* references = buffers.references;
    auto* vecmembers = buffers.vectorMembers;

    if (m_switch.isOn(collName)) {
      // Reconnect branches and collections
      m_datatree->SetBranchAddress(collName.c_str(), data);
      auto* colls = references;
      if (colls != nullptr) {
        for (size_t j = 0; j < colls->size(); ++j) {
          const auto brName = podio::root_utils::refBranch(collName, j);
          auto* l_branch = m_datatree->GetBranch(brName.c_str());
          l_branch->SetAddress(&(*colls)[j]);
        }
      }
      auto* colls_v = vecmembers;
      if (colls_v != nullptr) {
        int j = 0;
        for (auto& [dataType, add] : (*colls_v)) {
          const auto brName = podio::root_utils::vecBranch(collName, j);
          m_datatree->SetBranchAddress(brName.c_str(), add);
          ++j;
        }
      }
    }
    collBuffers->prepareForWrite();
  }
}

void PodioOutput::createBranches(const std::vector<std::pair<std::string, podio::CollectionBase*>>& collections) {
  for (const auto& [collName, collBuffers] : collections) {
    auto buffers = collBuffers->getBuffers();
    auto* data = buffers.data;
    auto* references = buffers.references;
    auto* vecmembers = buffers.vectorMembers;

    const std::string className     = collBuffers->getValueTypeName();
    const std::string collClassName = "vector<" + className + "Data>";

    int isOn = 0;
    if (m_switch.isOn(collName)) {
      isOn = 1;
      m_datatree->Branch(collName.c_str(), collClassName.c_str(), data);
      // Create branches for collections holding relations
      if (auto* refColls = references) {
        int j = 0;
        for (auto& c : (*refColls)) {
          const auto brName = podio::root_utils::refBranch(collName, j);
          m_datatree->Branch(brName.c_str(), c.get());
          ++j;
        }
      }
      // vector members
      if (auto* vminfo = vecmembers) {
        int j = 0;
        for (auto& [dataType, add] : (*vminfo)) {
          const std::string typeName = "vector<" + dataType + ">";
          const auto brName          = podio::root_utils::vecBranch(collName, j);
          m_datatree->Branch(brName.c_str(), typeName.c_str(), add);
          ++j;
        }
      }
    }

    const auto collID = m_podioDataSvc->getCollectionIDs()->collectionID(collName);
    const auto collType = collBuffers->getValueTypeName() + "Collection";
    m_collectionInfo.emplace_back(collID, std::move(collType), collBuffers->isSubsetCollection());

    debug() << isOn << " Registering collection " << collClassName << " " << collName.c_str() << " containing type "
            << className << endmsg;
    collBuffers->prepareForWrite();
  }
}

StatusCode PodioOutput::execute() {
  // for now assume identical content for every event
  // register for writing
  if (m_firstEvent) {
    m_collectionInfo.clear();
    createBranches(m_podioDataSvc->getCollections());
    createBranches(m_podioDataSvc->getReadCollections());
    m_metadatatree->Branch("CollectionTypeInfo", &m_collectionInfo);
  } else {
    resetBranches(m_podioDataSvc->getCollections());
    resetBranches(m_podioDataSvc->getReadCollections());
  }
  m_firstEvent = false;
  if (msgLevel(MSG::DEBUG)) {
    debug() << "Filling DataTree .." << endmsg;
  }
  m_datatree->Fill();
  m_evtMDtree->Fill();
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
  //// finalize trees and file //////////////////////////////
  debug() << "Finalizing trees and output file" << endmsg;
  m_file->cd();
  m_metadatatree->Branch("gaudiConfigOptions", &config_data);
  m_metadatatree->Branch("CollectionIDs", m_podioDataSvc->getCollectionIDs());
  m_metadatatree->Fill();
  m_colMDtree->Branch("colMD", "std::map<int,podio::GenericParameters>", m_podioDataSvc->getProvider().getColMetaDataMap() ) ;
  m_colMDtree->Fill();
  m_runMDtree->Branch("runMD", "std::map<int,podio::GenericParameters>", m_podioDataSvc->getProvider().getRunMetaDataMap() ) ;
  m_runMDtree->Fill();
  m_datatree->Write();
  m_file->Write();
  m_file->Close();
  info() << "Data written to: " << m_filename.value() << endmsg;
  if (!m_filenameRemote.value().empty()) {
    TFile::Cp(m_filename.value().c_str(), m_filenameRemote.value().c_str(), false);
    info() << " and copied to: " << m_filenameRemote.value() << endmsg;
  }
  return StatusCode::SUCCESS;
}
