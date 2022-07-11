#include "PileupHitMergeTool.h"

#include "podio/EventStore.h"

#include "edm4hep/CaloHitCollection.h"
#include "edm4hep/PositionedCaloHitCollection.h"
#include "edm4hep/PositionedTrackHitCollection.h"
#include "edm4hep/TrackHitCollection.h"

typedef PileupHitMergeTool<edm4hep::CaloHitCollection, edm4hep::PositionedCaloHitCollection> PileupCaloHitMergeTool;
typedef PileupHitMergeTool<edm4hep::TrackHitCollection, edm4hep::PositionedTrackHitCollection> PileupTrackHitMergeTool;
DECLARE_COMPONENT(PileupTrackHitMergeTool)
DECLARE_COMPONENT_WITH_ID(PileupTrackHitMergeTool, "PileupTrackHitMergeTool")
DECLARE_COMPONENT(PileupCaloHitMergeTool)
DECLARE_COMPONENT_WITH_ID(PileupCaloHitMergeTool, "PileupCaloHitMergeTool")

template <class Hits, class PositionedHits>
PileupHitMergeTool<Hits, PositionedHits>::PileupHitMergeTool(const std::string& aType, const std::string& aName,
                                                             const IInterface* aParent)
    : GaudiTool(aType, aName, aParent) {
  declareInterface<IEDMMergeTool>(this);
  declareProperty("signalHits", m_hitsSignal);
  declareProperty("signalPositionedHits", m_posHitsSignal);
  declareProperty("mergedHits", m_hitsMerged);
  declareProperty("mergedPositionedHits", m_posHitsMerged);
}

template <class Hits, class PositionedHits>
StatusCode PileupHitMergeTool<Hits, PositionedHits>::initialize() {
  return StatusCode::SUCCESS;
}

template <class Hits, class PositionedHits>
StatusCode PileupHitMergeTool<Hits, PositionedHits>::finalize() {
  return StatusCode::SUCCESS;
}

template <class Hits, class PositionedHits>
StatusCode PileupHitMergeTool<Hits, PositionedHits>::readPileupCollection(podio::EventStore& store) {
  // local pointers, to be filled by the event store
  const Hits* hitCollection;
  const PositionedHits* posHitCollection;

  // get collection address and store it in container
  bool hitCollectionPresent = store.get(m_pileupHitsBranchName, hitCollection);
  if (hitCollectionPresent) {
    m_hitCollections.push_back(hitCollection);
  } else {
    warning() << "No collection could be read from branch " << m_pileupHitsBranchName << endmsg;
    return StatusCode::FAILURE;
  }

  /// as above, for the positioned collection
  bool posHitCollectionPresent = store.get(m_pileupPosHitsBranchName, posHitCollection);
  if (posHitCollectionPresent) {
    m_posHitCollections.push_back(posHitCollection);
  } else {
    warning() << "No collection could be read from branch " << m_pileupPosHitsBranchName << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

template <class Hits, class PositionedHits>
StatusCode PileupHitMergeTool<Hits, PositionedHits>::readSignal() {
  // get collection from event sture
  auto collHitsSig = m_hitsSignal.get();
  auto collPosHitsSig = m_posHitsSignal.get();

  // store them in internal container
  m_hitCollections.push_back(collHitsSig);
  m_posHitCollections.push_back(collPosHitsSig);

  return StatusCode::SUCCESS;
}

template <class Hits, class PositionedHits>
StatusCode PileupHitMergeTool<Hits, PositionedHits>::mergeCollections() {

  // ownership given to data service at end of execute
  Hits* collHitsMerged = new Hits();
  PositionedHits* collPosHitsMerged = new PositionedHits();

  unsigned int collectionCounter = 0;
  for (auto hitColl : m_hitCollections) {
    // copy hits
    for (const auto elem : *hitColl) {

      auto clon = elem.clone();
      // add pileup vertex counter with an offset
      // i.e. for the signal event, 'bits' is just the trackID taken from geant
      // for the n-th pileup event, 'bits' is the trackID + n * offset
      // offset needs to be big enough to ensure uniqueness of trackID
      if (elem.bits() > m_trackIDCollectionOffset) {
        error() << "Event contains too many tracks to guarantee a unique trackID";
        error() << " The offset width or trackID field size needs to be adjusted!" << endmsg;
        return StatusCode::FAILURE;
      }

      clon.bits(clon.bits() + collectionCounter * m_trackIDCollectionOffset);
      collHitsMerged->push_back(clon);
    }
    ++collectionCounter;
  }
  for (auto posHitColl : m_posHitCollections) {
    // copy positioned hits
    for (const auto elem : *posHitColl) {
      collPosHitsMerged->push_back(elem.clone());
    }
  }

  m_hitsMerged.put(collHitsMerged);
  m_posHitsMerged.put(collPosHitsMerged);

  m_hitCollections.clear();
  m_posHitCollections.clear();
  return StatusCode::SUCCESS;
}
