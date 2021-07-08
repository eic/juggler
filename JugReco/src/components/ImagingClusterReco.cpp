/*
 *  Reconstruct the cluster/layer info for imaging calorimeter
 *  Logarithmic weighting is used to describe energy deposit in transverse direction
 *
 *  Author: Chao Peng (ANL), 06/02/2021
 */
#include "fmt/format.h"
#include <Eigen/Dense>
#include <algorithm>

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"

#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"
#include "JugBase/Utilities/Utils.hpp"

// Event Model related classes
#include "eicd/ImagingClusterCollection.h"
#include "eicd/ImagingLayerCollection.h"
#include "eicd/ImagingPixelCollection.h"

using namespace Gaudi::Units;
using namespace Eigen;

namespace Jug::Reco {

  class ImagingClusterReco : public GaudiAlgorithm {
  public:
    Gaudi::Property<double> m_sampFrac{this, "samplingFraction", 1.0};
    Gaudi::Property<int>    m_trackStopLayer{this, "trackStopLayer", 9};

    DataHandle<eic::ImagingPixelCollection> m_inputHitCollection{"inputHitCollection",
                                                                 Gaudi::DataHandle::Reader, this};
    DataHandle<eic::ImagingLayerCollection> m_outputLayerCollection{
        "outputLayerCollection", Gaudi::DataHandle::Writer, this};
    DataHandle<eic::ImagingClusterCollection> m_outputClusterCollection{
        "outputClusterCollection", Gaudi::DataHandle::Reader, this};

    ImagingClusterReco(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputLayerCollection", m_outputLayerCollection, "");
      declareProperty("outputClusterCollection", m_outputClusterCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }

      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& hits = *m_inputHitCollection.get();
      // output collections
      auto& layers   = *m_outputLayerCollection.createAndPut();
      auto& clusters = *m_outputClusterCollection.createAndPut();

      // Create a map of clusterID --> associated hits by looping over our clustered hits
      std::map<int, std::vector<eic::ConstImagingPixel>> cluster_map;
      for (const auto& hit : hits) {
        const size_t cid = hit.clusterID();
        if (!cluster_map.count(cid)) {
          cluster_map[cid] = {};
        }
        cluster_map[cid].push_back(hit);
      }

      for (const auto& [cid, clhits] : cluster_map) {
        // get cluster and associated layers
        auto cl        = reconstruct_cluster(clhits, cid);
        auto cl_layers = reconstruct_cluster_layers(clhits, cid);

        // reconstruct cluster direction
        const auto [cl_theta, cl_phi] = fit_track(cl_layers, m_trackStopLayer);
        cl.cl_theta(cl_theta);
        cl.cl_phi(cl_phi);

        // store layer and clusters on the datastore
        for (auto& layer : cl_layers) {
          layers.push_back(layer);
          // cl.addlayers(layer); // deprectated
        }
        clusters.push_back(cl);
      }

      // debug output
      for (const auto& cl : clusters) {
        debug() << fmt::format("Cluster {:d}: Edep = {:.3f} MeV, Dir = ({:.3f}, {:.3f}) deg",
                               cl.clusterID(), cl.edep()*1000., cl.cl_theta() / M_PI * 180.,
                               cl.cl_phi() / M_PI * 180.)
                << endmsg;
      }

      return StatusCode::SUCCESS;
    }

  private:
    template <typename T>
    static inline T pow2(const T& x)
    {
      return x * x;
    }

    std::vector<eic::ImagingLayer>
    reconstruct_cluster_layers(const std::vector<eic::ConstImagingPixel>& hits, const int cid) const
    {
      // using map to have hits sorted by layer
      std::map<int, std::vector<eic::ConstImagingPixel>> layer_map;
      for (const auto& hit : hits) {
        auto lid = hit.layerID();
        if (!layer_map.count(lid)) {
          layer_map[lid] = {};
        }
        layer_map[lid].push_back(hit);
      }

      // create layers
      std::vector<eic::ImagingLayer> cl_layers;
      for (const auto& [lid, lhits] : layer_map) {
        auto layer = reconstruct_layer(lhits, cid, lid);
        cl_layers.push_back(layer);
      }
      return cl_layers;
    }

    eic::ImagingLayer reconstruct_layer(const std::vector<eic::ConstImagingPixel>& hits,
                                        const int cid, const int lid) const
    {
      // use full members initialization here so it could catch changes in ecid
      eic::ImagingLayer layer{cid, lid, static_cast<int>(hits.size()), 0., 0., 0., 0., {}, {}};

      // mean position and total edep
      double mx   = 0.;
      double my   = 0.;
      double mz   = 0.;
      double edep = 0.;
      for (const auto& hit : hits) {
        mx += hit.x();
        my += hit.y();
        mz += hit.z();
        edep += hit.edep();
      }

      layer.x(mx / layer.nhits());
      layer.y(my / layer.nhits());
      layer.z(mz / layer.nhits());
      layer.edep(edep);

      double radius = 0.;
      for (const auto& hit : hits) {
        radius += std::sqrt(pow2(hit.x() - layer.x()) + pow2(hit.y() - layer.y()) +
                            pow2(hit.z() - layer.z()));
      }
      layer.radius(radius / layer.nhits());
      return layer;
    }

    eic::ImagingCluster reconstruct_cluster(const std::vector<eic::ConstImagingPixel>& hits,
                                            const int                                  cid) const
    {
      eic::ImagingCluster cluster;
      cluster.clusterID(cid);
      // eta, phi center, weighted by energy
      double meta = 0.;
      double mphi = 0.;
      double edep = 0.;
      double r    = 9999 * cm;
      for (const auto& hit : hits) {
        meta += hit.eta() * hit.edep();
        mphi += hit.phi() * hit.edep();
        edep += hit.edep();
        r = std::min(hit.r(), r);
      }
      const double eta   = meta / edep;
      const double phi   = mphi / edep;
      const double theta = 2. * std::atan(std::exp(-eta));
      cluster.nhits(hits.size());
      cluster.edep(edep);
      cluster.energy(edep / m_sampFrac); // simple energy reconstruction
      cluster.eta(eta);
      cluster.phi(phi);
      cluster.theta(theta);
      // project it to the inner surface (the shallowest hit)
      cluster.r(r);
      // cartesian coordinates
      ROOT::Math::Polar3DVectorD polar{r, theta, (phi > M_PI ? phi - M_PI : phi)};
      cluster.x(polar.x());
      cluster.y(polar.y());
      cluster.z(polar.z());

      // shower radius estimate (eta-phi plane)
      double radius = 0.;
      for (auto hit : cluster.hits()) {
        radius += std::sqrt(pow2(hit.eta() - cluster.eta()) + pow2(hit.phi() - cluster.phi()));
      }
      cluster.radius(radius / cluster.nhits());

      return cluster;
    }

    std::pair<double /* theta */, double /* phi */>
    fit_track(const std::vector<eic::ImagingLayer>& layers, const int stop_layer) const
    {
      int    nrows = 0;
      double mx    = 0.;
      double my    = 0.;
      double mz    = 0.;
      for (const auto& layer : layers) {
        if ((layer.layerID() <= stop_layer) && (layer.nhits() > 0)) {
          mx += layer.x();
          my += layer.y();
          mz += layer.z();
          nrows += 1;
        }
      }
      // cannot fit
      if (nrows < 2) {
        return {};
      }

      mx /= nrows;
      my /= nrows;
      mz /= nrows;
      // fill position data
      MatrixXd pos(nrows, 3);
      int      ir = 0;
      for (const auto& layer : layers) {
        if ((layer.layerID() <= stop_layer) && (layer.nhits() > 0)) {
          pos(ir, 0) = layer.x() - mx;
          pos(ir, 1) = layer.y() - my;
          pos(ir, 2) = layer.z() - mz;
          ir += 1;
        }
      }

      JacobiSVD<MatrixXd> svd(pos, ComputeThinU | ComputeThinV);
      // debug() << pos << endmsg;
      // debug() << svd.matrixV() << endmsg;
      const auto dir = svd.matrixV().col(0);
      // theta and phi
      return {std::acos(dir(2)), std::atan2(dir(1), dir(0))};
    }
  };

  DECLARE_COMPONENT(ImagingClusterReco)

} // namespace Jug::Reco

