#ifndef JugTrack_ProtoTrack_HH
#define JugTrack_ProtoTrack_HH

#include <cstddef>
#include <vector>

namespace Jug {

/// A proto track is a collection of hits identified by their indices.
using ProtoTrack = std::vector<size_t>;
/// Container of proto tracks. Each proto track is identified by its index.
using ProtoTrackContainer = std::vector<ProtoTrack>;

}  // namespace FW

#endif
