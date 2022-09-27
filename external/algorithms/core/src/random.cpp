#include <algorithms/random.h>

namespace algorithms {

void RandomSvc::init() {
  if (m_seed.hasValue()) {
    info() << "Custom random seed requested: " << m_seed << endmsg;
    m_gen = createBitGenerator(m_seed);
  }
}
void RandomSvc::init(const RandomSvc::GenFunc& gen) {
  init();
  info() << "Loading external generator function." << endmsg;
  m_gen = gen;
  if (m_seed.hasValue()) {
    warning() << "Custom random seed request ignored when using external generator function"
              << endmsg;
  }
}
RandomSvc::GenFunc RandomSvc::createBitGenerator(const size_t seed) {
  return [=](const size_t size) {
    static std::mutex m;
    static std::mt19937_64 gen{seed};
    std::lock_guard<std::mutex> lock{m};
    std::vector<CachedBitGenerator::value_type> ret{size};
    std::generate(ret.begin(), ret.end(), gen);
    return ret;
  };
}
} // namespace algorithms
