
#include <map>

#include <math.h>

#include "IRTAlgorithmServices.h"

#include <IRT/ParametricSurface.h>

// -------------------------------------------------------------------------------------

void Jug::Reco::IRTAlgorithmServices::configure_QE_lookup_table(const std::vector<std::pair<double, double>> &QE_vector, 
								    unsigned nbins)
{
  //const std::vector<std::pair<double, double>> &QE_vector = m_QE_input_data.value();
  //unsigned nbins = m_QEbins.value();

  // Well, could have probably just reordered the initial vector;
  std::map<double, double> QE_map;

  for(auto entry: QE_vector)
    QE_map[entry.first] = entry.second;

  // Sanity checks;
  if (QE_map.size() < 2 || nbins < 2) return;

  double from = (*QE_map.begin()).first, to = (*QE_map.rbegin()).first;
  // Will be "nbins+1" equidistant entries;
  double step = (to - from) / nbins;

  // Just in case somebody considers to call this method twice :-);
  m_QE_lookup_table.clear();

  for(auto entry: QE_map) {
    double e1 = entry.first, qe1 = entry.second;

    if (!m_QE_lookup_table.size())
      m_QE_lookup_table.push_back(std::make_pair(e1, qe1));
    else {
      const auto &prev = m_QE_lookup_table[m_QE_lookup_table.size()-1];

      double e0 = prev.first, qe0 = prev.second;
      double a = (qe1 - qe0) / (e1 - e0), b = qe0 - a*e0;
      // FIXME: check floating point accuracy when moving to a next point; do we actually 
      // care whether the overall number of bins will be "nbins+1" or more?;
      for(double e = e0+step; e<e1; e+=step)
	m_QE_lookup_table.push_back(std::make_pair(e, a*e + b));
    } //if
  } //for entry
} // Jug::Reco::IRTAlgorithmServices::configure_QE_lookup_table()

// -------------------------------------------------------------------------------------

bool Jug::Reco::IRTAlgorithmServices::QE_pass(double ev, double rand) const
{
  // Get the tabulated table reference; perform sanity checks;
  //const std::vector<std::pair<double, double>> &qe = u_quantumEfficiency.value();
  unsigned dim = m_QE_lookup_table.size(); if (dim < 2) return false;

  // Find a proper bin; no tricks, they are all equidistant;
  auto const &from = m_QE_lookup_table[0], &to = m_QE_lookup_table[dim-1];
  double emin = from.first, emax = to.first, step = (emax - emin) / (dim - 1);
  int ibin = (int)floor((ev - emin) / step);
  
  // Out of range check;
  if (ibin < 0 || ibin >= int(dim)) return false;

  // Get tabulated QE value, compare against the provided random variable;
  return (rand <= m_QE_lookup_table[ibin].second);
} // Jug::Reco::IRTAlgorithmServices::QE_pass()

// -------------------------------------------------------------------------------------

double Jug::Reco::IRTAlgorithmServices::GetDistance(const ParametricSurface *surface, 
							const eic::TrajectoryPoint *point) const
{
  if (!surface || !point) return 0.0;

  auto pt = TVector3(point->position.x, point->position.y, point->position.z);
  
  return surface->GetDistance(TVector3(point->position.x, point->position.y, point->position.z));
} // Jug::Reco::IRTAlgorithmServices::GetDistance()

// -------------------------------------------------------------------------------------

bool Jug::Reco::IRTAlgorithmServices::GetCrossing(const ParametricSurface *surface, 
						      const eic::TrajectoryPoint *point, 
						      TVector3 *crs) const
{
  if (!surface || !point) return false;
  
  return surface->GetCrossing(GetLocation(point), GetMomentum(point).Unit(), crs);
} // Jug::Reco::IRTAlgorithmServices::GetCrossing()

// -------------------------------------------------------------------------------------

TVector3 Jug::Reco::IRTAlgorithmServices::GetLocation(const eic::TrajectoryPoint *point) const
{
  return TVector3(point->position.x, point->position.y, point->position.z);
} // Jug::Reco::IRTAlgorithmServices::GetLocation()

// -------------------------------------------------------------------------------------

TVector3 Jug::Reco::IRTAlgorithmServices::GetMomentum(const eic::TrajectoryPoint *point) const
{
  double theta = point->p.theta, phi = point->p.phi; 

  return point->p.r*TVector3(sin(theta)*cos(phi), 
			     sin(theta)*sin(phi), 
			     cos(theta));
} // Jug::Reco::IRTAlgorithmServices::GetMomentum()

// -------------------------------------------------------------------------------------
