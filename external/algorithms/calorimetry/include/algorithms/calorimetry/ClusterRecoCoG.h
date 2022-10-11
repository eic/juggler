// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sylvester Joosten, Chao Peng, Whitney Armstrong

/*
 *  Reconstruct the cluster with Center of Gravity method
 *  Logarithmic weighting is used for mimicing energy deposit in transverse direction
 *
 *  Author: Sylvester Joosten (ANL), Chao Peng (ANL) 09/19/2022
 */

#include <algorithms/algorithm.h>
#include <algorithms/geo.h>

// Data types
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>
#include <edm4eic/ProtoClusterCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

namespace algorithms::calorimetry {

using ClusteringAlgorithm = Algorithm<
    Input<edm4eic::ProtoClusterCollection, std::optional<edm4hep::SimCalorimeterHitCollection>>,
    Output<edm4eic::ClusterCollection,
           std::optional<edm4eic::MCRecoClusterParticleAssociationCollection>>>;

/** Clustering with center of gravity method.
 *
 *  Reconstruct the cluster with Center of Gravity method
 *  Logarithmic weighting is used for mimicking energy deposit in transverse direction
 *
 * \ingroup reco
 */
class ClusterRecoCoG : public ClusteringAlgorithm {
public:
  using WeightFunc = std::function<double(double, double, double)>;

  // TODO: get rid of "Collection" in names
  ClusterRecoCoG(std::string_view name)
      : ClusteringAlgorithm{name,
                            {"inputProtoClusterCollection", "mcHits"},
                            {"outputClusterCollection", "outputAssociations"},
                            "Reconstruct a cluster with the Center of Gravity method. For "
                            "simulation results it optionally creates a Cluster <-> MCParticle "
                            "association provided both optional arguments are provided."} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  edm4eic::MutableCluster reconstruct(const edm4eic::ProtoCluster&) const;

  // TODO FIXME does the sampling fraction belong here or in the hit reconstruction?
  Property<double> m_sampFrac{this, "samplingFraction", 1.0, "Sampling fraction"};
  Property<double> m_logWeightBase{this, "logWeightBase", 3.6, "Weight base for log weighting"};
  Property<std::string> m_energyWeight{this, "energyWeight", "log", "Default hit weight method"};
  Property<std::string> m_moduleDimZName{this, "moduleDimZName", "", "z-dim name of the module"};
  // Constrain the cluster position eta to be within
  // the eta of the contributing hits. This is useful to avoid edge effects
  // for endcaps.
  Property<bool> m_enableEtaBounds{this, "enableEtaBounds", true, "Constrain cluster to hit eta?"};

  WeightFunc m_weightFunc;

  const GeoSvc& m_geo = GeoSvc::instance();
};
} // namespace algorithms::calorimetry

