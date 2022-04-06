// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Barak Schmookler

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

#include "JugBase/Utilities/Beam.h"
#include "JugBase/Utilities/Boost.h"

#include "Math/Vector4D.h"
using ROOT::Math::PxPyPzEVector;

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "eicd/MCRecoParticleAssociationCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/InclusiveKinematicsCollection.h"

namespace Jug::Reco {

class InclusiveKinematicsSigma : public GaudiAlgorithm {
private:
  DataHandle<edm4hep::MCParticleCollection> m_inputMCParticleCollection{
    "inputMCParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eicd::ReconstructedParticleCollection> m_inputParticleCollection{
    "inputReconstructedParticles",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eicd::MCRecoParticleAssociationCollection> m_inputParticleAssociation{
    "inputParticleAssociations",
    Gaudi::DataHandle::Reader,
    this};
  DataHandle<eicd::InclusiveKinematicsCollection> m_outputInclusiveKinematicsCollection{
    "outputInclusiveKinematics",
    Gaudi::DataHandle::Writer,
    this};

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton{0}, m_neutron{0}, m_electron{0};

public:
  InclusiveKinematicsSigma(const std::string& name, ISvcLocator* svcLoc)
      : GaudiAlgorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "MCParticles");
    declareProperty("inputReconstructedParticles", m_inputParticleCollection, "ReconstructedParticles");
    declareProperty("inputParticleAssociations", m_inputParticleAssociation, "MCRecoParticleAssociation");
    declareProperty("outputInclusiveKinematics", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsSigma");
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
    const auto ei_coll = Jug::Base::Beam::find_first_beam_electron(mcparts);
    if (ei_coll.size() == 0) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No beam electron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    const PxPyPzEVector ei(
      Jug::Base::Beam::round_beam_four_momentum(
        ei_coll[0].getMomentum(),
        m_electron,
        {-5.0, -10.0, -18.0},
        0.0)
      );

    // Get incoming hadron beam
    const auto pi_coll = Jug::Base::Beam::find_first_beam_hadron(mcparts);
    if (pi_coll.size() == 0) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No beam hadron found" << endmsg;
      }
      return StatusCode::SUCCESS;
    }
    const PxPyPzEVector pi(
      Jug::Base::Beam::round_beam_four_momentum(
        pi_coll[0].getMomentum(),
        pi_coll[0].getPDG() == 2212 ? m_proton : m_neutron,
        {41.0, 100.0, 275.0},
        m_crossingAngle)
      );

    // Get first scattered electron
    const auto ef_coll = Jug::Base::Beam::find_first_scattered_electron(mcparts);
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

    // Sums in colinear frame
    double pxsum = 0;
    double pysum = 0;
    double pzsum = 0;
    double Esum = 0;

    double pt_e = 0;
    double sigma_e = 0;

    // Get boost to colinear frame
    auto boost = Jug::Base::Boost::determine_boost(ei, pi);

    for(const auto& p: rcparts) {
      // Get the scattered electron index and angle
      if (p.id() == ef_rc_id) {
        // Lorentz vector in lab frame
        PxPyPzEVector e_lab(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
        // Boost to colinear frame
        PxPyPzEVector e_boosted = Jug::Base::Boost::apply_boost(boost, e_lab);

        pt_e = e_boosted.Pt();
        sigma_e = e_boosted.E() - e_boosted.Pz();

      // Sum over all particles other than scattered electron
      } else{
        // Lorentz vector in lab frame
        PxPyPzEVector hf_lab(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
        // Boost to colinear frame
        PxPyPzEVector hf_boosted = Jug::Base::Boost::apply_boost(boost, hf_lab);
        pxsum += hf_boosted.Px();
        pysum += hf_boosted.Py();
        pzsum += hf_boosted.Pz();
        Esum += hf_boosted.E();
      }
    }

    // DIS kinematics calculations
    auto sigma_h = Esum - pzsum;
    auto sigma_tot = sigma_e + sigma_h;

    if (sigma_h <= 0) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "No scattered electron found or sigma zero or negative" << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    // Calculate kinematic variables
    const auto y_sig = sigma_h / sigma_tot;
    const auto Q2_sig = (pt_e*pt_e) / (1. - y_sig);
    const auto x_sig = Q2_sig / (4.*ei.energy()*pi.energy()*y_sig);
    const auto nu_sig = Q2_sig / (2.*m_proton*x_sig);
    const auto W_sig = sqrt(m_proton*m_proton + 2*m_proton*nu_sig - Q2_sig);
    auto kin = out_kinematics.create(x_sig, Q2_sig, W_sig, y_sig, nu_sig);
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
DECLARE_COMPONENT(InclusiveKinematicsSigma)

} // namespace Jug::Reco
