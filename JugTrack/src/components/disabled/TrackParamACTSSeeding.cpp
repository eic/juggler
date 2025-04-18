#include <cmath>
#include <algorithm>
#include <unordered_map>

#include "Acts/ActsVersion.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinder.hpp"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include <k4FWCore/DataHandle.h>
#include <k4Interface/IGeoSvc.h>
#include "JugTrack/IActsGeoSvc.h"
#include "JugTrack/DD4hepBField.h"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include "DDRec/CellIDPositionConverter.h"

#include "edm4eic/TrackerHitCollection.h"
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
    class TrackParamACTSSeeding : public Gaudi::Algorithm {
    public:
        DataHandle<edm4eic::TrackerHitCollection>
        m_inputHitCollection { "inputHitCollection",
            Gaudi::DataHandle::Reader, this };
        DataHandle<ActsExamples::TrackParametersContainer>
        m_outputInitialTrackParameters {
            "outputInitialTrackParameters",
            Gaudi::DataHandle::Writer, this };

        SmartIF<IGeoSvc> m_geoSvc;
        SmartIF<IActsGeoSvc> m_actsGeoSvc;
        std::shared_ptr<const dd4hep::rec::CellIDPositionConverter> m_converter;

        // Acts::GeometryContext m_geoContext;
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
        class SpacePoint
        // : edm4eic::TrackerHitData
        {
        public:
            // Acts::Vector3 _dummy[16];
            Acts::Vector3 _position;
            Acts::Vector3 _positionError;
            int32_t _measurementIndex;
            const Acts::Surface *_surface;
            // Constructor to circumvent the fact that eic::TrackerHit
            // and associated classes are all non-polymorphic
            SpacePoint(const edm4eic::TrackerHit h,
                       const int32_t measurementIndex,
                       SmartIF<IActsGeoSvc> m_actsGeoSvc,
                       std::shared_ptr<const dd4hep::rec::CellIDPositionConverter> m_converter)
                : _measurementIndex(measurementIndex)
            {
                _position[0] = h.getPosition().x;
                _position[1] = h.getPosition().y;
                _position[2] = h.getPosition().z;
                _positionError[0] = h.getPositionError().xx;
                _positionError[1] = h.getPositionError().yy;
                _positionError[2] = h.getPositionError().zz;
                const auto volumeId =
                    m_converter->findContext(h.getCellID())->identifier;
                const auto its = m_actsGeoSvc->surfaceMap().find(volumeId);
                if (its == m_actsGeoSvc->surfaceMap().end()) {
                    _surface = nullptr;
                }
                else {
                    _surface = its->second;
                }
            }
            SpacePoint(const SpacePoint &sp)
                : _position(sp._position),
                  _positionError(sp._positionError),
                  _measurementIndex(sp._measurementIndex),
                  _surface(sp._surface)
            {
            }
            float x() const { return _position[0]; }
            float y() const { return _position[1]; }
            float z() const { return _position[2]; }
            float r() const { return std::hypot(x(), y()); }
            float varianceR() const
            {
                return (std::pow(x(), 2) * _positionError[0] +
                        std::pow(y(), 2) * _positionError[1]) /
                    (std::pow(x(), 2) + std::pow(y(), 2));
            }
            float varianceZ() const { return _positionError[2]; }
            constexpr uint32_t measurementIndex() const {
                return _measurementIndex; }
            bool isOnSurface() const {
                if (_surface == nullptr) {
                    return false;
                }
                return _surface->isOnSurface(
                     Acts::GeometryContext(), {x(), y(), z()},
                     {0, 0, 0});
            }
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
        using ProtoTrack = std::vector<ActsExamples::Index>;
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

            float cotThetaMax = std::sinh(4.01);
            float minPt = 100 * Acts::UnitConstants::MeV / cotThetaMax;
            float rMax = 440 * Acts::UnitConstants::mm;
            float zMin = -1500 * Acts::UnitConstants::mm;
            float zMax = 1700 * Acts::UnitConstants::mm;
            float deltaRMin = 0 * Acts::UnitConstants::mm;
            float deltaRMax = 600 * Acts::UnitConstants::mm;
            //
            float collisionRegionMin = -250 * Acts::UnitConstants::mm;
            float collisionRegionMax = 250 * Acts::UnitConstants::mm;
            float maxSeedsPerSpM = 2;
            float sigmaScattering = 5;
            float radLengthPerSeed = 0.1;
            float impactMax = 3 * Acts::UnitConstants::mm;
            //
            float bFieldInZ = 1.7 * Acts::UnitConstants::T;
            float beamPosX = 0 * Acts::UnitConstants::mm;
            float beamPosY = 0 * Acts::UnitConstants::mm;

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

            float deltaRMiddleMinSPRange = 10. * Acts::UnitConstants::mm;
            float deltaRMiddleMaxSPRange = 10. * Acts::UnitConstants::mm;

            // vector containing the map of z bins in the top and bottom layers
            std::vector<std::pair<int, int> > zBinNeighborsTop;
            std::vector<std::pair<int, int> > zBinNeighborsBottom;
        } m_cfg;
        Acts::SpacePointGridConfig m_gridCfg;
        Acts::SpacePointGridOptions m_gridOpt;
        Acts::SeedFinderConfig<SpacePoint> m_finderCfg;
        Acts::SeedFinderOptions m_finderOpt;
        /// The track parameters covariance (assumed to be the same
        /// for all estimated track parameters for the moment)
        Acts::BoundSquareMatrix m_covariance =
            Acts::BoundSquareMatrix::Zero();

    public:
        TrackParamACTSSeeding(const std::string &name,
                              ISvcLocator* svcLoc)
            : Gaudi::Algorithm(name, svcLoc)
        {
            declareProperty("inputHitCollection",
                            m_inputHitCollection, "");
            declareProperty("outputInitialTrackParameters",
                            m_outputInitialTrackParameters, "");
        }

        StatusCode initialize() override;

        StatusCode execute(const EventContext&) const override;
    };


    StatusCode TrackParamACTSSeeding::initialize()
    {
        if (Gaudi::Algorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }

        m_geoSvc = service("GeoSvc");
        if (m_geoSvc == nullptr) {
            error() << "Unable to locate Geometry Service. " << endmsg;
            return StatusCode::FAILURE;
        }
        m_actsGeoSvc = service("ActsGeoSvc");
        if (m_actsGeoSvc == nullptr) {
            error() << "Unable to locate ACTS Geometry Service. " << endmsg;
            return StatusCode::FAILURE;
        }              
        m_converter = std::make_shared<const dd4hep::rec::CellIDPositionConverter>(*(m_geoSvc->getDetector()));

        m_BField = std::dynamic_pointer_cast<
            const Jug::BField::DD4hepBField>(
                m_actsGeoSvc->getFieldProvider());
        m_fieldContext = Jug::BField::BFieldVariant(m_BField);

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
        m_finderCfg.impactMax = m_cfg.impactMax;

        m_finderOpt.bFieldInZ = m_cfg.bFieldInZ;
        m_finderOpt.beamPos =
            Acts::Vector2(m_cfg.beamPosX, m_cfg.beamPosY);

        m_finderCfg =
          m_finderCfg.toInternalUnits().calculateDerivedQuantities();
        m_finderOpt =
          m_finderOpt.toInternalUnits().calculateDerivedQuantities(
            m_finderCfg);

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

    StatusCode TrackParamACTSSeeding::execute(const EventContext&) const
    {
        const edm4eic::TrackerHitCollection *hits =
            m_inputHitCollection.get();
        // Create output collections
        auto initTrackParameters =
            m_outputInitialTrackParameters.createAndPut();

        static SeedContainer seeds;
        static Acts::SeedFinder<SpacePoint>::SeedingState state;

        // Sadly, eic::TrackerHit and eic::TrackerHitData are
	// non-polymorphic
        std::vector<SpacePoint> spacePoint;
        std::vector<const SpacePoint *> spacePointPtrs;

        std::shared_ptr<const Acts::TrackingGeometry>
            trackingGeometry = m_actsGeoSvc->trackingGeometry();

#ifdef USE_LOCAL_COORD
        // Currently broken, possibly geometry issues
#error Do not use, broken
        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__ << ": "
                    << sourceLinks->size() << ' '
                    << hits->size() << endmsg;
        }
        auto its = sourceLinks->begin();
        auto itm = measurements->begin();
        for (; its != sourceLinks->end() &&
                 itm != measurements->end();
             its++, itm++) {
            const Acts::Surface *surface =
                trackingGeometry->findSurface(its->get().geometryId());
            if (surface != nullptr) {
                Acts::Vector3 v =
                    surface->localToGlobal(
                        Acts::GeometryContext(),
                        {std::get<Acts::Measurement<
                         Acts::BoundIndices, 2>>(*itm).
                         parameters()[0],
                         std::get<Acts::Measurement<
                         Acts::BoundIndices, 2>>(*itm).
                         parameters()[1]},
                        {0, 0, 0});
                if (msgLevel(MSG::DEBUG)) {
                    debug() << __FILE__ << ':' << __LINE__ << ": "
                            << its - sourceLinks->begin() << ' '
                        // << itm - measurements->begin() << ' '
                            << v[0] << ' ' << v[1] << ' ' << v[2]
                            << endmsg;
                }
                spacePoint.push_back(
                    SpacePoint(
                        edm4eic::TrackerHit(
                            static_cast<uint64_t>(spacePoint.size() + 1),
                            edm4eic::Vector3f(v[0], v[1], v[2]),
                            // Dummy uncertainties
                            edm4eic::CovDiag3f(25.0e-6 / 3.0,
                                               25.0e-6 / 3.0, 0.0),
                            0.0, 10.0, 0.05, 0.0),
                        static_cast<int32_t>(spacePoint.size())));
                spacePointPtrs.push_back(&spacePoint.back());
        }
#else // USE_LOCAL_COORD
        for(const auto &h : *hits) {
            spacePoint.push_back(SpacePoint(
                h, static_cast<int32_t>(spacePoint.size() + 1),
                m_actsGeoSvc, m_converter));
            spacePointPtrs.push_back(&spacePoint.back());
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
                        << ' ' << spacePointPtrs.back()
                        << ' ' << spacePointPtrs.back()->measurementIndex()
                        << ' ' << spacePointPtrs.back()->isOnSurface()
                        << endmsg;
            }
        }
#endif // USE_LOCAL_COORD
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

        // extent used to store r range for middle spacepoint
        Acts::Extent rRangeSPExtent;

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

        const float bFieldInZSave = m_gridOpt.bFieldInZ;
        const float minPtSave = m_gridCfg.minPt;
        m_gridCfg.minPt = 400 * Acts::UnitConstants::MeV;
        m_gridOpt.bFieldInZ =
            (m_gridCfg.minPt / Acts::UnitConstants::MeV) /
            (150.0 * (1.0 + FLT_EPSILON) *
             (m_gridCfg.rMax / Acts::UnitConstants::mm)) *
            1000.0 * Acts::UnitConstants::T;
        if (msgLevel(MSG::DEBUG)) {
            debug() << "createGrid() with temporary minPt = "
                    << m_gridCfg.minPt / Acts::UnitConstants::MeV
                    << " MeV, bFieldInZ = "
                    << m_gridOpt.bFieldInZ / (1000 * Acts::UnitConstants::T)
                    << " kT" << endmsg;
        }
        auto grid =
            Acts::SpacePointGridCreator::createGrid<SpacePoint>(
                m_gridCfg, m_gridOpt);
        m_gridOpt.bFieldInZ = bFieldInZSave;
        m_gridCfg.minPt = minPtSave;
        if (msgLevel(MSG::DEBUG)) {
            debug() << "phiBins = "
                    << grid->axes().front()->getNBins()
                    << ", zBins = "
                    << grid->axes().back()->getNBins() << endmsg;
        }

        auto spacePointsGrouping =
            Acts::BinnedSPGroup<SpacePoint>(
                spacePointPtrs.begin(), spacePointPtrs.end(),
                extractGlobalQuantities, bottomBinFinder,
                topBinFinder, std::move(grid),
                rRangeSPExtent,
                m_finderCfg, m_finderOpt);
        auto finder = Acts::SeedFinder<SpacePoint>(m_finderCfg);

        if (msgLevel(MSG::DEBUG)) {
            debug() << __FILE__ << ':' << __LINE__
                    << ": spacePointsGrouping.size() = "
                    << spacePointsGrouping.size() << endmsg;
        }

        const Acts::Range1D<float> rMiddleSPRange(
            std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
            m_cfg.deltaRMiddleMinSPRange,
            std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2 -
            m_cfg.deltaRMiddleMaxSPRange);

        // Run the seeding
        seeds.clear();

        for (const auto [bottom, middle, top] : spacePointsGrouping) {
            finder.createSeedsForGroup(
                m_finderOpt, state, spacePointsGrouping.grid(),
                std::back_inserter(seeds), bottom, middle, top, rMiddleSPRange
            );
        }

        if (msgLevel(MSG::DEBUG)) {
            debug() << "seeds.size() = " << seeds.size() << endmsg;
        }

        ActsExamples::TrackParametersContainer trackParameters;
        ProtoTrackContainer tracks;
        trackParameters.reserve(seeds.size());
        tracks.reserve(seeds.size());

        std::shared_ptr<const Acts::MagneticFieldProvider>
            magneticField = m_actsGeoSvc->getFieldProvider();

        // if (msgLevel(MSG::DEBUG)) { debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg; }
        auto bCache = magneticField->makeCache(m_fieldContext);

        std::unordered_map<size_t, bool> spTaken;

        for (size_t iseed = 0; iseed < seeds.size(); iseed++) {
            const auto &seed = seeds[iseed];
            // Get the bottom space point and its reference surface
            const auto bottomSP = seed.sp().front();
            auto hitIdx = bottomSP->measurementIndex();
            // const Acts::Surface *surface = nullptr;
            const Acts::Surface *surface = bottomSP->_surface;
            if (surface == nullptr) {
                if (msgLevel(MSG::DEBUG)) {
                    debug() << "hit " << hitIdx << " ("
                            << bottomSP->x() << ", " << bottomSP->y()
                            << ", " << bottomSP->z()
                            << ") lost its surface" << endmsg;
                }
            }
            if (std::find_if(spacePoint.begin(), spacePoint.end(),
                             [&surface](const SpacePoint sp) {
                                 return surface == sp._surface;
                             }) == spacePoint.end()) {
                if (msgLevel(MSG::DEBUG)) {
                    debug() << "hit " << hitIdx
                            << " has a surface that was never "
                            << "associated with a hit" << endmsg;
                }
                surface = nullptr;
            }
            if (surface == nullptr && msgLevel(MSG::DEBUG)) {
                debug() << "hit " << hitIdx
                        << " is not found in the tracking gemetry"
                        << endmsg;
                continue;
            }

            // Get the magnetic field at the bottom space point
            auto fieldRes = magneticField->getField(
                {bottomSP->x(), bottomSP->y(), bottomSP->z()},
                bCache);
            // if (msgLevel(MSG::DEBUG)) { debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg; }
            // Estimate the track parameters from seed
            auto optParams = Acts::estimateTrackParamsFromSeed(
                Acts::GeometryContext(),
                seed.sp().begin(), seed.sp().end(),
                *surface, *fieldRes, m_cfg.bFieldMin);
            // if (msgLevel(MSG::DEBUG)) { debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg; }
            if (not optParams.has_value()) {
                debug() << "Estimation of track parameters for seed "
                        << iseed << " failed." << endmsg;
                continue;
            }
            else if (!(spTaken[seed.sp()[0]->measurementIndex()] ||
                       spTaken[seed.sp()[1]->measurementIndex()] ||
                       spTaken[seed.sp()[2]->measurementIndex()])) {
                const auto& params = optParams.value();
                initTrackParameters->emplace_back(
                    surface->getSharedPtr(), params,
                    m_covariance, Acts::ParticleHypothesis::pion());
                // Activate/deactivate for unique seed filtering
#if 0
                spTaken[seed.sp()[0]->measurementIndex()] = true;
                spTaken[seed.sp()[1]->measurementIndex()] = true;
                spTaken[seed.sp()[2]->measurementIndex()] = true;
#endif
            }
            // if (msgLevel(MSG::DEBUG)) { debug() << __FILE__ << ':' << __LINE__ << ": " << endmsg; }
        }

        return StatusCode::SUCCESS;
    }
    
    DECLARE_COMPONENT(TrackParamACTSSeeding)
} // namespace Jug::reco
