
#include <vector>

#include <TVector3.h>

#include "DDRec/CellIDPositionConverter.h"

#ifndef _IRT_ALGORITHM_SERVICES_
#define _IRT_ALGORITHM_SERVICES_

#include "eicd/TrajectoryPoint.h"
class ParametricSurface;

namespace Jug::PID {

  // Wanto to split off the part, which knows nothing about Gaudi & Co; somehow this is 
  // less trivial than expected (need to pass the initial parameters separately);
  class IRTAlgorithmServices {
  public:
    IRTAlgorithmServices( void ) {};
    ~IRTAlgorithmServices( void ) {};

  protected:
    std::vector<std::pair<double, double>> m_QE_lookup_table;

    std::vector<std::pair<double, double>> ApplyFineBinning(const std::vector<std::pair<double, double>> &input, 
							    unsigned nbins);
    bool GetFinelyBinnedTableEntry(const std::vector<std::pair<double, double>> &table, 
				   double argument, double *value) const;

    // Keep the same method name as in Chao's PhotoMultiplierDigi.cpp; however implement 
    // a fast equidistant array lookup instead;
    bool QE_pass(double ev, double rand) const;

    // FIXME: this crap needs to be optimized;
    double GetDistance(const ParametricSurface *surface, const eic::TrajectoryPoint *point) const;
    bool GetCrossing(const ParametricSurface *surface, const eic::TrajectoryPoint *point, 
		     TVector3 *crs) const;
    TVector3 GetLocation(const eic::TrajectoryPoint *point) const;
    TVector3 GetMomentum(const eic::TrajectoryPoint *point) const;

    void findPrimitive(
        const std::string typeName, const dd4hep::Solid sol,
        dd4hep::Solid &prim, dd4hep::Position &pos, const TGeoMatrix *matx=nullptr
        ) const;
  };

} // namespace Jug::PID

#endif
