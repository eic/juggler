#include "Gaudi/Algorithm.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Producer.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/PhysicalConstants.h"
#include <algorithm>
#include <cmath>

#include "JugBase/IParticleSvc.h"
#include "JugBase/DataHandle.h"

#include "JugReco/Utilities/Beam.h"
#include "JugReco/Utilities/Boost.h"

#include "Math/Vector4D.h"
using ROOT::Math::PxPyPzEVector;

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/MCRecoParticleAssociationCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Reco {

class InclusiveKinematicsDA : public GaudiAlgorithm {
private:
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticleCollection{
    "MCParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eicd::ReconstructedParticleCollection> m_inputParticleCollection{
    "ReconstructedParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eicd::MCRecoParticleAssociationCollection> m_inputParticleAssociation{
    "MCRecoParticleAssociation",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eicd::InclusiveKinematicsCollection> m_outputInclusiveKinematicsCollection{
    "InclusiveKinematicsDA",
    Gaudi::DataHandle::Writer,
    this};

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton{0}, m_neutron{0}, m_electron{0};

public:
  InclusiveKinematicsDA(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "MCParticles");
    declareProperty("inputParticles", m_inputParticleCollection, "ReconstructedParticles");
    declareProperty("inputAssociations", m_inputParticleAssociation, "MCRecoParticleAssociation");
    declareProperty("outputData", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsDA");
  }

  StatusCode initialize() override {
    if (GaudiAlgorithm::initialize().isFailure())
      return StatusCode::FAILURE;

    m_pidSvc = service("ParticleSvc");
    if (!m_pidSvc) {
      error() << "Unable to locate Particle Service. "
              << "Make sure you have ParticleSvc in the configuration."
              << endmsg;
      return StatusCode::FAILURE;
    }
    m_proton = m_pidSvc->particle(2212).mass;
    m_neutron = m_pidSvc->particle(2112).mass;
    m_electron = m_pidSvc->particle(11).mass;

    return StatusCode::SUCCESS;
  }

  StatusCode execute() override {
    // input collections
    const auto& mcparts = *(m_inputMCParticleCollection.get());
    const auto& rcparts = *(m_inputParticleCollection.get());
    const auto& rcassoc = *(m_inputParticleAssociation.get());
    // output collection
    auto& out_kinematics = *(m_outputInclusiveKinematicsCollection.createAndPut());

    // Get incoming electron beam
    const auto ei_coll = Jug::Reco::Beam::find_first_beam_electron(mcparts);
    if (ei_coll.size() == 0) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No beam electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    const PxPyPzEVector ei(
      Jug::Reco::Beam::round_beam_four_momentum(
        ei_coll[0].getMomentum(),
        m_electron,
        {-5.0, -10.0, -18.0},
        0.0)
      );

    // Get incoming hadron beam
    const auto pi_coll = Jug::Reco::Beam::find_first_beam_hadron(mcparts);
    if (pi_coll.size() == 0) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No beam hadron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    const PxPyPzEVector pi(
      Jug::Reco::Beam::round_beam_four_momentum(
        pi_coll[0].getMomentum(),
        pi_coll[0].getPDG() == 2212 ? m_proton : m_neutron,
        {41.0, 100.0, 275.0},
        m_crossingAngle)
      );

    // Get first scattered electron
    const auto ef_coll = Jug::Reco::Beam::find_first_scattered_electron(mcparts);
    if (ef_coll.size() == 0) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No truth scattered electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    // Associate first scattered electron with reconstructed electrons
    //const auto ef_assoc = std::find_if(
    //  rcassoc.begin(),
    //  rcassoc.end(),
    //  [&ef_coll](const auto& a){ return a.getSimID() == ef_coll[0].id(); });
    auto ef_assoc = rcassoc.begin();
    for (; ef_assoc != rcassoc.end(); ++ef_assoc) {
      if (ef_assoc->getSimID() == ef_coll[0].id()) {
        break;
      }
    }
    if (!(ef_assoc != rcassoc.end())) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "Truth scattered electron not in reconstructed particles" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    const auto ef_rc{ef_assoc->getRec()};
    const auto ef_rc_id{ef_rc.id()};

    // Loop over reconstructed particles to get all outgoing particles
    // ----------------------------------------------------------------- 
    // Right now, everything is taken from Reconstructed particles branches.
    // 
    // This means the tracking detector is used for charged particles to caculate the momentum,
    // and the magnitude of this momentum plus the true PID to calculate the energy.
    // No requirement is made that these particles produce a hit in any other detector
    // 
    // Using the Reconstructed particles branches also means that the reconstruction for neutrals is done using the
    // calorimeter(s) information for the energy and angles, and then using this energy and the true PID to get the
    // magnitude of the momentum.
    // -----------------------------------------------------------------

    //Sums in colinear frame
    double pxsum = 0;
    double pysum = 0;
    double pzsum = 0;
    double Esum = 0;
    double theta_e = 0;

    // Get boost to colinear frame
    auto boost = Jug::Reco::Boost::determine_boost(ei, pi);

    for(const auto& p : rcparts) {
      // Get the scattered electron index and angle
      if (p.id() == ef_rc_id) {
        // Lorentz vector in lab frame
        PxPyPzEVector e_lab(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
        // Boost to colinear frame
        PxPyPzEVector e_boosted = Jug::Reco::Boost::apply_boost(boost, e_lab);

        theta_e = e_boosted.Theta();

      // Sum over all particles other than scattered electron
      } else {
        // Lorentz vector in lab frame
        PxPyPzEVector hf_lab(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
        // Boost to colinear frame
        PxPyPzEVector hf_boosted = Jug::Reco::Boost::apply_boost(boost, hf_lab);

        pxsum += hf_boosted.Px();
        pysum += hf_boosted.Py();
        pzsum += hf_boosted.Pz();
        Esum += hf_boosted.E();
      }
    }

    // DIS kinematics calculations
    auto sigma_h = Esum - pzsum;
    auto ptsum = sqrt(pxsum*pxsum + pysum*pysum);
    auto theta_h = 2.*atan(sigma_h/ptsum);

    // If no scattered electron was found
    if (sigma_h <= 0) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No scattered electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    // Calculate kinematic variables
    const auto y_da = tan(theta_h/2.) / ( tan(theta_e/2.) + tan(theta_h/2.) );
    const auto Q2_da = 4.*ei.energy()*ei.energy() * ( 1. / tan(theta_e/2.) ) * ( 1. / (tan(theta_e/2.) + tan(theta_h/2.)) );
    const auto x_da = Q2_da / (4.*ei.energy()*pi.energy()*y_da);
    const auto nu_da = Q2_da / (2.*m_proton*x_da);
    const auto W_da = sqrt(m_proton*m_proton + 2*m_proton*nu_da - Q2_da);
    auto kin = out_kinematics.create(x_da, Q2_da, W_da, y_da, nu_da);
    kin.setScat(ef_rc);

    // Debugging output
    if (msgLevel(MSG::DEBUG)) {
      debug() << "pi = " << pi << endmsg;
      debug() << "ei = " << ei << endmsg;
      debug() << "x,Q2,W,y,nu = "
              << kin.getX() << ","
              << kin.getQ2() << ","
              << kin.getW() << ","
              << kin.getY() << ","
              << kin.getNu()
              << endmsg;
    }

    return StatusCode::SUCCESS;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
DECLARE_COMPONENT(InclusiveKinematicsDA)

} // namespace Jug::Reco
