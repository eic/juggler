#ifndef JutTrack_ProtoVertex_H
#define JutTrack_ProtoVertex_H

#include "JugTrack/Index.hpp"

#include <vector>

// Copied from ACTS Examples
//
namespace Jug {

/// A proto vertex is a collection of tracks identified by their indices.
using ProtoVertex = std::vector<Index>;
/// Container of proto vertices. Each proto vertex is identified by its index.
using ProtoVertexContainer = std::vector<ProtoVertex>;

}  // namespace ActsExamples

#endif
