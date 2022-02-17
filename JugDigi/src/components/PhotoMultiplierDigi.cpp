/*  General PhotoMultiplier Digitization
 *
 *  Apply the given quantum efficiency for photon detection
 *  Converts the number of detected photons to signal amplitude
 *
 *  Author: Chao Peng (ANL)
 *  Date: 10/02/2020
 */

#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <cmath>

#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/PhysicalConstants.h"

#include "JugBase/DataHandle.h"
#include "JugBase/UniqueID.h"

// Event Model related classes
#include "eicd/RawPMTHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"


using namespace Gaudi::Units;

namespace Jug::Digi {

/** PhotoMultiplierDigi.
 *
 * \ingroup digi
 */
class PhotoMultiplierDigi : public GaudiAlgorithm, AlgorithmIDMixin<>
{
public:
    DataHandle<edm4hep::TrackerHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::RawPMTHitCollection>
        m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
    Gaudi::Property<std::vector<std::pair<double, double>>>
        u_quantumEfficiency{this, "quantumEfficiency", {{2.6*eV, 0.3}, {7.0*eV, 0.3}}};
    Gaudi::Property<double> m_hitTimeWindow{this, "hitTimeWindow", 20.0*ns};
    Gaudi::Property<double> m_timeStep{this, "timeStep", 0.0625*ns};
    Gaudi::Property<double> m_speMean{this, "speMean", 80.0};
    Gaudi::Property<double> m_speError{this, "speError", 16.0};
    Gaudi::Property<double> m_pedMean{this, "pedMean", 200.0};
    Gaudi::Property<double> m_pedError{this, "pedError", 3.0};
    Rndm::Numbers m_rngUni, m_rngNorm;

    // constructor
    PhotoMultiplierDigi(const std::string& name, ISvcLocator* svcLoc)
        : GaudiAlgorithm(name, svcLoc)
        , AlgorithmIDMixin(name, info())
    {
        declareProperty("inputHitCollection", m_inputHitCollection,"");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
        if (GaudiAlgorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }

        auto randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        auto sc1 = m_rngUni.initialize(randSvc, Rndm::Flat(0., 1.));
        auto sc2 = m_rngNorm.initialize(randSvc, Rndm::Gauss(0., 1.));
        if (!sc1.isSuccess() || !sc2.isSuccess()) {
            error() << "Cannot initialize random generator!" << endmsg;
            return StatusCode::FAILURE;
        }

        qe_init();

        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // input collection
        const auto &sim = *m_inputHitCollection.get();
        // Create output collections
        auto &raw = *m_outputHitCollection.createAndPut();

        struct HitData { int npe; double signal; double time; };
        std::unordered_map<long long, std::vector<HitData>> hit_groups;
        // collect the photon hit in the same cell
        // calculate signal
        for(const auto& ahit : sim) {
            // quantum efficiency
            if (!qe_pass(ahit.energyDeposit(), m_rngUni())) {
                continue;
            }
            // cell id, time, signal amplitude
            long long id = ahit.cellID();
            double time = ahit.truth().time;
            double amp = m_speMean + m_rngNorm()*m_speError;

            // group hits
            auto it = hit_groups.find(id);
            if (it != hit_groups.end()) {
                size_t i = 0;
                for (auto git = it->second.begin(); git != it->second.end(); ++git, ++i) {
                    if (std::abs(time - git->time) <= (m_hitTimeWindow/ns)) {
                        git->npe += 1;
                        git->signal += amp;
                        break;
                    }
                }
                // no hits group found
                if (i >= it->second.size()) {
                    it->second.emplace_back(HitData{1, amp + m_pedMean + m_pedError*m_rngNorm(), time});
                }
            } else {
                hit_groups[id] = {HitData{1, amp + m_pedMean + m_pedError*m_rngNorm(), time}};
            }
        }

        // build hit
        int ID = 0;
        for (auto &it : hit_groups) {
            for (auto &data : it.second) {
                eic::RawPMTHit hit{
                  {ID++, algorithmID()},
                  it.first,
                  static_cast<uint32_t>(data.signal), 
                  static_cast<uint32_t>(data.time/(m_timeStep/ns))};
                raw.push_back(hit);
            }
        }

        return StatusCode::SUCCESS;
    }

private:
    void qe_init()
    {
        auto &qeff = u_quantumEfficiency.value();

        // sort quantum efficiency data first
        std::sort(qeff.begin(), qeff.end(),
            [] (const std::pair<double, double> &v1, const std::pair<double, double> &v2) {
                return v1.first < v2.first;
            });

        // sanity checks
        if (qeff.empty()) {
            qeff = {{2.6*eV, 0.3}, {7.0*eV, 0.3}};
            warning() << "Invalid quantum efficiency data provided, using default values: " << qeff << endmsg;
        }
        if (qeff.front().first > 3.0*eV) {
            warning() << "Quantum efficiency data start from " << qeff.front().first/eV
                      << " eV, maybe you are using wrong units?" << endmsg;
        }
        if (qeff.back().first < 6.0*eV) {
            warning() << "Quantum efficiency data end at " << qeff.back().first/eV
                      << " eV, maybe you are using wrong units?" << endmsg;
        }
    }

    // helper function for linear interpolation
    // Comp return is defined as: equal, 0;  greater, > 0; less, < 0
    template<class RndmIter, typename T, class Compare>
    RndmIter interval_search(RndmIter beg, RndmIter end, const T &val, Compare comp) const
    {
        // special cases
        auto dist = std::distance(beg, end);
        if ((dist < 2) || (comp(*beg, val) > 0) || (comp(*std::prev(end), val) < 0)) {
            return end;
        }
        auto mid = std::next(beg, dist / 2);

        while (mid != end) {
            if (comp(*mid, val) == 0) {
                return mid;
            } else if (comp(*mid, val) > 0) {
                end = mid;
            } else {
                beg = std::next(mid);
            }
            mid = std::next(beg, std::distance(beg, end)/2);
        }

        if (mid == end || comp(*mid, val) > 0) {
            return std::prev(mid);
        }
        return mid;
    }

    bool qe_pass(double ev, double rand) const
    {
        auto &qeff = u_quantumEfficiency.value();
        auto it = interval_search(qeff.begin(), qeff.end(), ev,
                    [] (const std::pair<double, double> &vals, double val) {
                        return vals.first - val;
                    });

        if (it == qeff.end()) {
            // info() << ev/eV << " eV is out of QE data range, assuming 0% efficiency" << endmsg;
            return false;
        }

        double prob = it->second;
        auto itn = std::next(it);
        if (itn != qeff.end() && (itn->first - it->first != 0)) {
            prob = (it->second*(itn->first - ev) + itn->second*(ev - it->first)) / (itn->first - it->first);
        }

        // info() << ev/eV << " eV, QE: "  << prob*100. << "%" << endmsg;
        return rand <= prob;
    }
};

DECLARE_COMPONENT(PhotoMultiplierDigi)

} // namespace Jug::Digi

