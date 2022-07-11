#include "PileupHitMergeTool.h"

#include "podio/EventStore.h"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

typedef PileupHitMergeTool<edm4hep::SimCalorimeterHitCollection> PileupCaloHitMergeTool;
typedef PileupHitMergeTool<edm4hep::SimTrackerHitCollection> PileupTrackHitMergeTool;
DECLARE_COMPONENT(PileupTrackHitMergeTool)
DECLARE_COMPONENT(PileupCaloHitMergeTool)

template <class Hits>
PileupHitMergeTool<Hits>::PileupHitMergeTool(const std::string& aType, const std::string& aName,
                                             const IInterface* aParent)
    : GaudiTool(aType, aName, aParent) {
  declareProperty("signalHits", m_hitsSignal);
  declareProperty("mergedHits", m_hitsMerged);
}

template <class Hits>
StatusCode PileupHitMergeTool<Hits>::initialize() {
  return StatusCode::SUCCESS;
}

template <class Hits>
StatusCode PileupHitMergeTool<Hits>::finalize() {
  return StatusCode::SUCCESS;
}

template <class Hits>
StatusCode PileupHitMergeTool<Hits>::readPileupCollection(podio::EventStore& store) {
  // local pointers, to be filled by the event store
  const Hits* hitCollection;

  // get collection address and store it in container
  bool hitCollectionPresent = store.get(m_pileupHitsBranchName, hitCollection);
  if (hitCollectionPresent) {
    m_hitCollections.push_back(hitCollection);
  } else {
    warning() << "No collection could be read from branch " << m_pileupHitsBranchName << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

template <class Hits>
StatusCode PileupHitMergeTool<Hits>::readSignal() {
  // get collection from event sture
  auto collHitsSig = m_hitsSignal.get();

  // store them in internal container
  m_hitCollections.push_back(collHitsSig);

  return StatusCode::SUCCESS;
}

template <class Hits>
StatusCode PileupHitMergeTool<Hits>::mergeCollections() {

  // ownership given to data service at end of execute
  Hits* collHitsMerged = new Hits();

  unsigned int collectionCounter = 0;
  for (auto hitColl : m_hitCollections) {
    // copy hits
    for (const auto elem : *hitColl) {

      auto clon = elem.clone();
      // add pileup vertex counter with an offset
      // i.e. for the signal event, 'bits' is just the trackID taken from geant
      // for the n-th pileup event, 'bits' is the trackID + n * offset
      // offset needs to be big enough to ensure uniqueness of trackID
/*
      if (elem.bits() > m_trackIDCollectionOffset) {
        error() << "Event contains too many tracks to guarantee a unique trackID";
        error() << " The offset width or trackID field size needs to be adjusted!" << endmsg;
        return StatusCode::FAILURE;
      }

      clon.bits(clon.bits() + collectionCounter * m_trackIDCollectionOffset);
*/
      collHitsMerged->push_back(clon);
    }
    ++collectionCounter;
  }

  m_hitsMerged.put(collHitsMerged);
  m_hitCollections.clear();
  return StatusCode::SUCCESS;
}
