#include "PileupParticlesMergeTool.h"

#include "podio/EventStore.h"

#include "edm4hep/GenVertexCollection.h"
#include "edm4hep/MCParticleCollection.h"

DECLARE_COMPONENT(PileupParticlesMergeTool)

PileupParticlesMergeTool::PileupParticlesMergeTool(const std::string& aType, const std::string& aName,
                                                   const IInterface* aParent)
    : GaudiTool(aType, aName, aParent) {
  declareInterface<IEDMMergeTool>(this);
  declareProperty("signalGenVertices", m_vertIn);
  declareProperty("signalGenParticles", m_partIn);
  declareProperty("mergedGenParticles", m_partOut);
  declareProperty("mergedGenVertices", m_vertOut);
}

StatusCode PileupParticlesMergeTool::initialize() { return StatusCode::SUCCESS; }

StatusCode PileupParticlesMergeTool::finalize() { return StatusCode::SUCCESS; }

StatusCode PileupParticlesMergeTool::readPileupCollection(podio::EventStore& store) {
  const edm4hep::MCParticleCollection* mcParticleCollection;
  const edm4hep::GenVertexCollection* genVertexCollection;

  bool mcParticleCollectionPresent = store.get(m_pileupGenParticlesBranchName, mcParticleCollection);
  debug() << "size of collection about to be read " << mcParticleCollection->size() << endmsg;
  if (mcParticleCollectionPresent) {
    m_MCParticleCollections.push_back(mcParticleCollection);
  } else {
    warning() << "No collection could be read from branch " << m_pileupGenParticlesBranchName << endmsg;
  }

  bool genVertexCollectionPresent = store.get(m_pileupGenVerticesBranchName, genVertexCollection);
  if (genVertexCollectionPresent) {
    m_GenVertexCollections.push_back(genVertexCollection);
  } else {
    warning() << "No collection could be read from branch " << m_pileupGenVerticesBranchName << endmsg;
  }

  return StatusCode::SUCCESS;
}

StatusCode PileupParticlesMergeTool::readSignal() {
  auto collVSig = m_vertIn.get();
  auto collPSig = m_partIn.get();

  m_MCParticleCollections.push_back(collPSig);
  m_GenVertexCollections.push_back(collVSig);

  return StatusCode::SUCCESS;
}

StatusCode PileupParticlesMergeTool::mergeCollections() {

  debug() << "merge " << m_GenVertexCollections.size() << "  collections ..." << endmsg;

  // ownership given to data service at end of execute
  edm4hep::MCParticleCollection* collPOut = new edm4hep::MCParticleCollection();
  edm4hep::GenVertexCollection* collVOut = new edm4hep::GenVertexCollection();

  // need to keep track of the accumulated length
  // of collections already merged
  unsigned int collectionOffset = 0;

  // merge vertices

  for (auto genVertexColl : m_GenVertexCollections) {
    debug() << " traverse collection of size " << genVertexColl->size() << endmsg;
    for (const auto elem : *genVertexColl) {
      collVOut->push_back(elem.clone());
    }
  }
  // merge particles, keeping the references up to date
  for (unsigned int collCounter = 0; collCounter < m_MCParticleCollections.size(); ++collCounter) {
    auto mcPartColl = m_MCParticleCollections[collCounter];
    for (const auto elem : *mcPartColl) {
      auto newPart = edm4hep::MCParticle(elem.core());
      if (elem.startVertex().isAvailable()) {
        // update reference: find new index in merged collection
        // add offset since the signal particles are copied in new collection first
        unsigned int newIdStart = elem.startVertex().getObjectID().index + collectionOffset;
        // update startVertex
        newPart.startVertex((*collVOut)[newIdStart]);
      }
      if (elem.endVertex().isAvailable()) {
        // update reference: find new index in merged collection
        // add offset since the signal particles are copied in new collection first
        unsigned int newIdEnd = elem.endVertex().getObjectID().index + collectionOffset;
        // update startVertex
        newPart.endVertex((*collVOut)[newIdEnd]);
      }
      collPOut->push_back(newPart);
    }
    collectionOffset += m_GenVertexCollections[collCounter]->size();
  }

  m_vertOut.put(collVOut);
  m_partOut.put(collPOut);

  m_MCParticleCollections.clear();
  m_GenVertexCollections.clear();
  return StatusCode::SUCCESS;
}
