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

// Event Model related classes
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/ClusterCollection.h"

using namespace Gaudi::Units;

namespace Jug::Reco {

  class EcalTungstenSamplingCluster : public GaudiAlgorithm {
  public:
    Gaudi::Property<double>                   m_minModuleEdep{this, "minModuleEdep", 0.5 * MeV};
    DataHandle<eic::CalorimeterHitCollection> m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader,
                                                                   this};
    DataHandle<eic::ClusterCollection> m_outputClusterCollection{"outputClusterCollection", Gaudi::DataHandle::Writer,
                                                                 this};

    /// Pointer to the geometry service
    SmartIF<IGeoSvc> m_geoSvc;

    // ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    EcalTungstenSamplingCluster(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
      declareProperty("inputHitCollection", m_inputHitCollection, "");
      declareProperty("outputClusterCollection", m_outputClusterCollection, "");
    }

    StatusCode initialize() override
    {
      if (GaudiAlgorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
      }
      m_geoSvc = service("GeoSvc");
      if (!m_geoSvc) {
        error() << "Unable to locate Geometry Service. "
                << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
        return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
      // input collections
      const auto& hits = *m_inputHitCollection.get();
      // Create output collections
      auto& clusters = *m_outputClusterCollection.createAndPut();

      // Loop over hits
      std::vector<bool> visits(hits.size(), false);
      double ref_pos_z = 0.0;
      for (size_t i = 0; i < hits.size(); i++)
      {
	// Ignore hits that already visited
	if (visits[i]) {
		continue;
	}
	// Select above minimum energy deposit
      	if (hits[i].energy() > m_minModuleEdep/GeV) {
          // Reference z position of hit
	  ref_pos_z = hits[i].position().z;
	  std::cout << "Before ref_pos_z: " << ref_pos_z << std::endl;
	  // Call function to add up energy deposit in the same layer based on z position
	  group_by_layer(i, ref_pos_z, hits, visits, clusters);
	}
      }
      return StatusCode::SUCCESS;
    }
  private:
    // Grouping hits by layers
    void group_by_layer(int index, double ref_pos_z, const eic::CalorimeterHitCollection &hits, std::vector<bool> &visits, eic::ClusterCollection &clusters) const
    {
	visits[index] = true;
	auto tot_edep = hits[index].energy();
	double pos_x = hits[index].position().x;
	double pos_y = hits[index].position().y;
	double pos_z = hits[index].position().z;
	double temp = ref_pos_z;
	std::cout << "After ref_pos_z: " << temp << std::endl;
	std::cout << "Reprint pos_z: " << pos_z << std::endl;
	// Loop over hits to find hits on the same layer
	for (size_t i = 0; i < hits.size(); i++) {
		if(visits[i]) {
			continue;
		}
		// Add up energy deposit based on the same z position and above energy threshold
		if((double)hits[i].position().z == ref_pos_z && hits[i].energy() > m_minModuleEdep/GeV) {
			tot_edep += hits[i].energy();
			visits[i] = true;
		}
	}
	// Save info as a cluster
	// TODO: position x and y determination
	clusters.push_back(eic::Cluster{tot_edep,{pos_x,pos_y,pos_z}, {0,0,0,0,0,0}, 0, 0});
	return;
    }
  };
  DECLARE_COMPONENT(EcalTungstenSamplingCluster)
} // namespace Jug::Reco
