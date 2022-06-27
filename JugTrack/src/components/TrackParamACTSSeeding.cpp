#include <cmath>
#include <algorithm>
#include <unordered_map>

#include "Acts/ActsVersion.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/BField/DD4hepBField.h"
#include "JugTrack/IndexSourceLink.hpp"
#include "JugTrack/Measurement.hpp"
#include "JugTrack/Track.hpp"

#include "eicd/TrackerHitCollection.h"

#include "Math/Vector3D.h"


  ///// (Reconstructed) track parameters e.g. close to the vertex.
  //using TrackParameters = Acts::CurvilinearTrackParameters;

  ///// Container of reconstructed track states for multiple tracks.
  //using TrackParametersContainer = std::vector<TrackParameters>;

  ///// MultiTrajectory definition
  //using Trajectory = Acts::MultiTrajectory<SourceLink>;

  ///// Container for the truth fitting/finding track(s)
  //using TrajectoryContainer = std::vector<SimMultiTrajectory>;

namespace Jug::Reco {

    /** Initial Track parameters from MC truth.
     *
     *  TrackParmetersContainer
     *  \ingroup tracking
     */
    class TrackParamACTSSeeding : public GaudiAlgorithm {
    public:
        DataHandle<IndexSourceLinkContainer>
        m_inputSourceLinks { "inputSourceLinks",
            Gaudi::DataHandle::Reader, this };
        DataHandle<MeasurementContainer>
        m_inputMeasurements { "inputMeasurements",
            Gaudi::DataHandle::Reader, this};
        DataHandle<eicd::TrackerHitCollection>
        m_inputHitCollection { "inputHitCollection",
            Gaudi::DataHandle::Reader, this };
        DataHandle<TrackParametersContainer>
        m_outputInitialTrackParameters {
            "outputInitialTrackParameters",
            Gaudi::DataHandle::Writer, this };

        SmartIF<IGeoSvc> m_geoSvc;
        Acts::GeometryContext m_geoContext;
        std::shared_ptr<const Jug::BField::DD4hepBField> m_BField =
            nullptr;
        Acts::MagneticFieldContext m_fieldContext;

        /// Index type to reference elements in a container.
        ///
        /// We do not expect to have more than 2^32 elements in any
        /// given container so a fixed sized integer type is
        /// sufficient.
        //using Index = eic::Index;

        /// Space point representation of eic::TrackerHitData suitable
        /// for ACTS track seeding.
        class SpacePoint : eicd::TrackerHitData {
        public:
            int32_t _measurementIndex;
            // Constructor to circumvent the fact that eic::TrackerHit
            // and associated classes are all non-polymorphic
            SpacePoint(const eicd::TrackerHit h,
                       const int32_t measurementIndex)
                : _measurementIndex(measurementIndex)
            {
                position = h.getPosition();
                positionError = h.getPositionError();
            }
            constexpr float x() const { return position.x; }
            constexpr float y() const { return position.y; }
            constexpr float z() const { return position.z; }
            constexpr float r() const { return std::hypot(x(), y()); }
            constexpr float varianceR() const
            {
                return (std::pow(x(), 2) * positionError.xx +
                        std::pow(y(), 2) * positionError.yy) /
                    (std::pow(x(), 2) + std::pow(y(), 2));
            }
            constexpr float varianceZ() const { return positionError.zz; }
            constexpr uint32_t measurementIndex() const {
                return _measurementIndex; }
        };

        static bool spCompare(SpacePoint r, SpacePoint s)
        {
            return
                std::hypot(r.x(), r.y(), r.z()) <
                std::hypot(s.x(), s.y(), s.z());
        }

        /// Container of sim seed
        using SeedContainer = std::vector<Acts::Seed<SpacePoint>>;

        /// A proto track is a collection of hits identified by their
        /// indices.
        using ProtoTrack = std::vector<Index>;
        /// Container of proto tracks. Each proto track is identified by
        /// its index.
        using ProtoTrackContainer = std::vector<ProtoTrack>;

        using SpacePointContainer = std::vector<SpacePoint>;

        struct Config {
            /// Input space point collections.
            ///
            /// We allow multiple space point collections to allow
            /// different parts of the detector to use different
            /// algorithms for space point construction, e.g.
            /// single-hit space points for pixel-like detectors or
            /// double-hit space points for strip-like detectors.
            std::vector<std::string> inputSpacePoints;
            /// Output track seed collection.
            std::string outputSeeds;
            /// Output proto track collection.
            std::string outputProtoTracks;

            float bFieldInZ = 3 * Acts::UnitConstants::T;
            float minPt = 500 * Acts::UnitConstants::MeV;
            float rMax = 440 * Acts::UnitConstants::mm;
            float zMin = -1500 * Acts::UnitConstants::mm;
            float zMax = 1700 * Acts::UnitConstants::mm;
            float deltaRMin = 1 * Acts::UnitConstants::mm;
            float deltaRMax = 600 * Acts::UnitConstants::mm;
            float cotThetaMax = sinh(4.01);
            //
            float collisionRegionMin = -250 * Acts::UnitConstants::mm;
            float collisionRegionMax = 250 * Acts::UnitConstants::mm;
            float maxSeedsPerSpM = 6;
            float sigmaScattering = 5;
            float radLengthPerSeed = 0.1;
            float beamPosX = 0 * Acts::UnitConstants::mm;
            float beamPosY = 0 * Acts::UnitConstants::mm;
            float impactMax = 3 * Acts::UnitConstants::mm;

            /// The minimum magnetic field to trigger the track
            /// parameters estimation
            double bFieldMin = 0.1 * Acts::UnitConstants::T;

            /// Constant term of the loc0 resolution.
            double sigmaLoc0 = 25 * Acts::UnitConstants::um;
            /// Constant term of the loc1 resolution.
            double sigmaLoc1 = 100 * Acts::UnitConstants::um;
            /// Phi angular resolution.
            double sigmaPhi = 0.02 * Acts::UnitConstants::degree;
            /// Theta angular resolution.
            double sigmaTheta = 0.02 * Acts::UnitConstants::degree;
            /// q/p resolution.
            double sigmaQOverP = 0.1 / Acts::UnitConstants::GeV;
            /// Time resolution.
            double sigmaT0 = 1400 * Acts::UnitConstants::s;

	    int numPhiNeighbors = 1;

            // vector containing the map of z bins in the top and bottom layers
            std::vector<std::pair<int, int> > zBinNeighborsTop;
            std::vector<std::pair<int, int> > zBinNeighborsBottom;
        } m_cfg;
        Acts::SpacePointGridConfig m_gridCfg;
        Acts::SeedfinderConfig<SpacePoint> m_finderCfg;
        /// The track parameters covariance (assumed to be the same
        /// for all estimated track parameters for the moment)
        Acts::BoundSymMatrix m_covariance =
            Acts::BoundSymMatrix::Zero();

    public:
        TrackParamACTSSeeding(const std::string &name,
                              ISvcLocator* svcLoc)
            : GaudiAlgorithm(name, svcLoc) {
            declareProperty("inputSourceLinks",
                            m_inputSourceLinks, "");
            declareProperty("inputMeasurements",
                            m_inputMeasurements, "");
            declareProperty("inputHitCollection",
                            m_inputHitCollection, "");
            declareProperty("outputInitialTrackParameters",
                            m_outputInitialTrackParameters, "");
        }

        StatusCode initialize() override;

        void
        findSeed(SeedContainer &seeds,
                 const eicd::TrackerHitCollection *hits,
                 const IndexSourceLinkContainer *sourceLinks,
                 const MeasurementContainer *measurements,
                 Acts::Seedfinder<SpacePoint>::State &state);

        StatusCode execute() override;
    };


    StatusCode TrackParamACTSSeeding::initialize()
    {
        if (GaudiAlgorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }

        m_geoSvc = service("GeoSvc");
        if (m_geoSvc == nullptr) {
            error() << "Unable to locate Geometry Service. " << endmsg;
            return StatusCode::FAILURE;
        }

        m_BField = std::dynamic_pointer_cast<
            const Jug::BField::DD4hepBField>(
                m_geoSvc->getFieldProvider());
        m_fieldContext = Jug::BField::BFieldVariant(m_BField);

        m_gridCfg.bFieldInZ = m_cfg.bFieldInZ;
        m_gridCfg.minPt = m_cfg.minPt;
        m_gridCfg.rMax = m_cfg.rMax;
        m_gridCfg.zMax = m_cfg.zMax;
        m_gridCfg.zMin = m_cfg.zMin;
        m_gridCfg.cotThetaMax = m_cfg.cotThetaMax;

        m_gridCfg.deltaRMax = m_cfg.deltaRMax;

        // Construct seed filter
        Acts::SeedFilterConfig filterCfg;
        filterCfg.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
        m_finderCfg.seedFilter =
            std::make_unique<Acts::SeedFilter<SpacePoint>>(
                Acts::SeedFilter<SpacePoint>(filterCfg));

        m_finderCfg.rMax = m_cfg.rMax;
        m_finderCfg.deltaRMin = m_cfg.deltaRMin;
        m_finderCfg.deltaRMax = m_cfg.deltaRMax;
        m_finderCfg.collisionRegionMin = m_cfg.collisionRegionMin;
        m_finderCfg.collisionRegionMax = m_cfg.collisionRegionMax;
        m_finderCfg.zMin = m_cfg.zMin;
        m_finderCfg.zMax = m_cfg.zMax;
        m_finderCfg.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
        m_finderCfg.cotThetaMax = m_cfg.cotThetaMax;
        m_finderCfg.sigmaScattering = m_cfg.sigmaScattering;
        m_finderCfg.radLengthPerSeed = m_cfg.radLengthPerSeed;
        m_finderCfg.minPt = m_cfg.minPt;
        m_finderCfg.bFieldInZ = m_cfg.bFieldInZ;
        m_finderCfg.beamPos =
            Acts::Vector2(m_cfg.beamPosX, m_cfg.beamPosY);
        m_finderCfg.impactMax = m_cfg.impactMax;

        // Set up the track parameters covariance (the same for all
        // tracks)
        m_covariance(Acts::eBoundLoc0, Acts::eBoundLoc0) =
            std::pow(m_cfg.sigmaLoc0, 2);
        m_covariance(Acts::eBoundLoc1, Acts::eBoundLoc1) =
            std::pow(m_cfg.sigmaLoc1, 2);
        m_covariance(Acts::eBoundPhi, Acts::eBoundPhi) =
            std::pow(m_cfg.sigmaPhi, 2);
        m_covariance(Acts::eBoundTheta, Acts::eBoundTheta) =
            std::pow(m_cfg.sigmaTheta, 2);
        m_covariance(Acts::eBoundQOverP, Acts::eBoundQOverP) =
            std::pow(m_cfg.sigmaQOverP, 2);
        m_covariance(Acts::eBoundTime, Acts::eBoundTime) =
            std::pow(m_cfg.sigmaT0, 2);

        return StatusCode::SUCCESS;
    }

    void TrackParamACTSSeeding::
    findSeed(SeedContainer &seeds,
             const eicd::TrackerHitCollection *hits,
             const IndexSourceLinkContainer *sourceLinks,
             const MeasurementContainer *measurements,
             Acts::Seedfinder<SpacePoint>::State &state)
    {
        // Sadly, eic::TrackerHit and eic::TrackerHitData are
	// non-polymorphic
        std::vector<SpacePoint> spacePoint;
        std::vector<const SpacePoint *> spacePointPtrs;
        // extent used to store r range for middle spacepoint
        Acts::Extent rRangeSPExtent;

        std::shared_ptr<const Acts::TrackingGeometry>
            trackingGeometry = m_geoSvc->trackingGeometry();

        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__ << ": "
                    << sourceLinks->size() << ' '
                    << measurements->size() << ' '
                    << hits->size() << endmsg;
        }
        auto its = sourceLinks->begin();
        auto itm = measurements->begin();
        for (; its != sourceLinks->end() &&
                 itm != measurements->end();
             its++, itm++) {
            const Acts::Surface *surface = trackingGeometry->findSurface(its->get().geometryId());
            if (surface != nullptr) {
                Acts::Vector3 v = surface->localToGlobal(m_geoContext, {std::get<Acts::Measurement<Acts::BoundIndices, 2>>(*itm).parameters()[0], std::get<Acts::Measurement<Acts::BoundIndices, 2>>(*itm).parameters()[1]}, {0, 0, 0});
                if (msgLevel(MSG::DEBUG)) {
                    debug() << __FILE__ << ':' << __LINE__ << ": "
                            << its - sourceLinks->begin() << ' '
                        // << itm - measurements->begin() << ' '
                            << v[0] << ' ' << v[1] << ' ' << v[2]
                            << endmsg;
                }
#ifdef USE_LOCAL_COORD
                spacePoint.push_back(
                    SpacePoint(
                        eicd::TrackerHit(
                            static_cast<uint64_t>(spacePoint.size()),
                            eicd::Vector3f(v[0], v[1], v[2]),
                            eicd::CovDiag3f(25.0e-6 / 3.0,
                                            25.0e-6 / 3.0, 0.0),
                            0.0, 10.0, 0.05, 0.0),
                        static_cast<int32_t>(spacePoint.size())));
                spacePointPtrs.push_back(&spacePoint.back());
                rRangeSPExtent.check({ spacePoint.back().x(),
                                       spacePoint.back().y(),
                                       spacePoint.back().z() });
#endif // USE_LOCAL_COORD
            }
        }

        for(const auto &h : *hits) {
            if (msgLevel(MSG::DEBUG)) {
                debug() << __FILE__ << ':' << __LINE__ << ": "
                        << ' ' << h.getPosition().x
                        << ' ' << h.getPosition().y
                        << ' ' << h.getPosition().z
                        << ' ' << h.getPositionError().xx
                        << ' ' << h.getPositionError().yy
                        << ' ' << h.getPositionError().zz
                        << ' ' << h.getTime()
                        << ' ' << h.getTimeError()
                        << ' ' << h.getEdep()
                        << ' ' << h.getEdepError()
                        << endmsg;
            }
#ifndef USE_LOCAL_COORD
            spacePoint.push_back(SpacePoint(h, static_cast<int32_t>(spacePoint.size())));
            spacePointPtrs.push_back(&spacePoint.back());
            rRangeSPExtent.check({ spacePoint.back().x(),
                                   spacePoint.back().y(),
                                   spacePoint.back().z() });
#endif // USE_LOCAL_COORD
        }
        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg;
        }

        auto extractGlobalQuantities =
            [=](const SpacePoint& sp, float, float, float) ->
            std::pair<Acts::Vector3, Acts::Vector2> {
                Acts::Vector3 position { sp.x(), sp.y(), sp.z() };
                Acts::Vector2 covariance {
                    sp.varianceR(), sp.varianceZ() };
                return std::make_pair(position, covariance);
        };
        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg;
        }

        auto bottomBinFinder =
            std::make_shared<Acts::BinFinder<SpacePoint>>(
                Acts::BinFinder<SpacePoint>(m_cfg.zBinNeighborsBottom,
                                            m_cfg.numPhiNeighbors));
        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg;
        }
        auto topBinFinder =
            std::make_shared<Acts::BinFinder<SpacePoint>>(
                Acts::BinFinder<SpacePoint>(m_cfg.zBinNeighborsTop,
                                            m_cfg.numPhiNeighbors));
        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg;
        }

        auto grid =
            Acts::SpacePointGridCreator::createGrid<SpacePoint>(
                m_gridCfg);
        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg;
        }

        auto spacePointsGrouping =
            Acts::BinnedSPGroup<SpacePoint>(
                spacePointPtrs.begin(), spacePointPtrs.end(),
                extractGlobalQuantities, bottomBinFinder,
                topBinFinder, std::move(grid), m_finderCfg);
        auto finder = Acts::Seedfinder<SpacePoint>(m_finderCfg);

        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__
                    << ": spacePointsGrouping.size() = "
                    << spacePointsGrouping.size() << endmsg;
        }
#if 0
        topBinFinder.get();
        bottomBinFinder.get();
#endif
        // Run the seeding
        seeds.clear();

        auto group = spacePointsGrouping.begin();
        auto groupEnd = spacePointsGrouping.end();
#if 1
        for (; !(group == groupEnd); ++group) {
            finder.createSeedsForGroup(
                state, std::back_inserter(seeds),
                group.bottom(), group.middle(), group.top(),
                rRangeSPExtent);
        }
#endif

        if (msgLevel(MSG::DEBUG)) {
            debug() << "seeds.size() = " << seeds.size() << endmsg;
        }
    }

    StatusCode TrackParamACTSSeeding::execute()
    {
        const eicd::TrackerHitCollection *hits =
            m_inputHitCollection.get();
        const IndexSourceLinkContainer *sourceLinks =
            m_inputSourceLinks.get();
        const MeasurementContainer *measurements =
            m_inputMeasurements.get();
        // Create output collections
        auto initTrackParameters =
            m_outputInitialTrackParameters.createAndPut();

        static thread_local SeedContainer seeds;
        static thread_local Acts::Seedfinder<SpacePoint>::State state;

        findSeed(seeds, hits, sourceLinks, measurements, state);

        TrackParametersContainer trackParameters;
        ProtoTrackContainer tracks;
        trackParameters.reserve(seeds.size());
        tracks.reserve(seeds.size());

        std::shared_ptr<const Acts::TrackingGeometry>
            trackingGeometry = m_geoSvc->trackingGeometry();
        std::shared_ptr<const Acts::MagneticFieldProvider>
            magneticField = m_geoSvc->getFieldProvider();

        if (msgLevel(MSG::DEBUG)) { debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg; }
        auto bCache = magneticField->makeCache(m_fieldContext);

        std::unordered_map<size_t, bool> spTaken;

        for (size_t iseed = 0; iseed < seeds.size(); iseed++) {
            const auto &seed = seeds[iseed];
            // Get the bottom space point and its reference surface
            const auto bottomSP = seed.sp().front();
            auto hitIdx = bottomSP->measurementIndex();
            // if (msgLevel(MSG::DEBUG)) {
            //     // debug() << __FILE__ << ':' << __LINE__ << ": iseed = " << iseed << ", seed.sp().size() = " << seed.sp().size() << ", hitIdx = " << hitIdx << endmsg;
            //     for (auto i : seed.sp()) {
            //         // debug() << __FILE__ << ':' << __LINE__ << ": " << i->measurementIndex() << ", " << i->x() << ", " << i->y() << ", " << i->z() << endmsg;
            //     }
            // }
            // Guard against any memory access issues
            hitIdx = std::min(hitIdx, static_cast<uint32_t>(
                sourceLinks->size() - 1));
            const Acts::Surface *surface = nullptr;
            for (auto &s : *sourceLinks) {
                surface = trackingGeometry->findSurface(s.get().geometryId());
                if (surface != nullptr &&
                    surface->isOnSurface(
                        m_geoContext,
                        {bottomSP->x(), bottomSP->y(), bottomSP->z()},
                        {0, 0, 0})) {
                    break;
                }
            }
            if (surface == nullptr && msgLevel(MSG::DEBUG)) {
                debug() << "hit " << hitIdx
                        << " is not found in the tracking gemetry"
                        << endmsg;
                continue;
            }

            if (msgLevel(MSG::DEBUG)) {
                debug() << __FILE__ << ':' << __LINE__
                        << ": iseed = " << iseed << ", "
                        << surface->type() << ", "
                        << surface->center(m_geoContext).x() << ", "
                        << surface->center(m_geoContext).y() << ", "
                        << surface->center(m_geoContext).z()
                        << endmsg;
            }

            // Get the magnetic field at the bottom space point
            auto fieldRes = magneticField->getField(
                {bottomSP->x(), bottomSP->y(), bottomSP->z()},
                bCache);
            if (msgLevel(MSG::DEBUG)) { debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg; }
            // Estimate the track parameters from seed
            auto optParams = Acts::estimateTrackParamsFromSeed(
                m_geoContext, seed.sp().begin(), seed.sp().end(),
                *surface, *fieldRes, m_cfg.bFieldMin);
            if (msgLevel(MSG::DEBUG)) { debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg; }
            if (not optParams.has_value()) {
                debug() << "Estimation of track parameters for seed "
                        << iseed << " failed." << endmsg;
                continue;
            }
            else if (!(spTaken[seed.sp()[0]->measurementIndex()] ||
                       spTaken[seed.sp()[1]->measurementIndex()] ||
                       spTaken[seed.sp()[2]->measurementIndex()])) {
                const auto& params = optParams.value();
                const double charge =
                    std::copysign(1, params[Acts::eBoundQOverP]);
                initTrackParameters->emplace_back(
                    surface->getSharedPtr(), params, charge,
                    m_covariance);
                // Activate/deactivate for unique seed filtering
#if 0
                spTaken[seed.sp()[0]->measurementIndex()] = true;
                spTaken[seed.sp()[1]->measurementIndex()] = true;
                spTaken[seed.sp()[2]->measurementIndex()] = true;
#endif
            }
            if (msgLevel(MSG::DEBUG)) { debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg; }
        }

        return StatusCode::SUCCESS;
    }
    
    DECLARE_COMPONENT(TrackParamACTSSeeding)
} // namespace Jug::reco
