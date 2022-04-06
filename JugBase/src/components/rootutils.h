// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Wouter Deconinck

#ifndef PODIO_ROOT_UTILS_H
#define PODIO_ROOT_UTILS_H

#include "podio/podioVersion.h"
#include "podio/CollectionBase.h"
#include "podio/CollectionBranches.h"

#include "TBranch.h"
#include "TClass.h"

#include <vector>
#include <string>

namespace podio::root_utils {
// Workaround slow branch retrieval for 6.22/06 performance degradation
// see: https://root-forum.cern.ch/t/serious-degradation-of-i-o-performance-from-6-20-04-to-6-22-06/43584/10
template<class Tree>
TBranch* getBranch(Tree* chain, const char* name) {
  return static_cast<TBranch*>(chain->GetListOfBranches()->FindObject(name));
}

inline std::string refBranch(const std::string& name, size_t index) {
  return name + "#" + std::to_string(index);
}

inline std::string vecBranch(const std::string& name, size_t index) {
  return name + "_" + std::to_string(index);
}


inline void setCollectionAddresses(podio::CollectionBase* collection, const CollectionBranches& branches) {
  auto buffers = collection->getBuffers();
  auto data = buffers.data;
  auto references = buffers.references;
  auto vecmembers = buffers.vectorMembers;

  if (data) {
    branches.data->SetAddress(data);
  }

  if (references) {
    for (size_t i = 0; i < references->size(); ++i) {
      branches.refs[i]->SetAddress(&(*references)[i]);
    }
  }

  if (vecmembers) {
    for (size_t i = 0; i < vecmembers->size(); ++i) {
      branches.vecs[i]->SetAddress((*vecmembers)[i].second);
    }
  }
}

}

#endif
