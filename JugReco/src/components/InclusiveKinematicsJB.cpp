// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck

#include "Gaudi/Algorithm.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/RndmGenerators.h"
#include <algorithm>
#include <cmath>

#include "JugBase/IParticleSvc.h"
#include <k4FWCore/DataHandle.h>

#include "JugBase/Utilities/Beam.h"
#include "JugBase/Utilities/Boost.h"

#include "Math/Vector4D.h"
using ROOT::Math::PxPyPzEVector;

// Event Model related classes
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/InclusiveKinematicsCollection.h"

namespace Jug::Reco {

class InclusiveKinematicsJB : public Gaudi::Algorithm {
private:
  mutable DataHandle<edm4hep::MCParticleCollection> m_inputMCParticleCollection{
    "inputMCParticles",
    Gaudi::DataHandle::Reader,
    this};
  mutable DataHandle<edm4eic::ReconstructedParticleCollection> m_inputParticleCollection{
    "inputReconstructedParticles",
    Gaudi::DataHandle::Reader,
    this};
  mutable DataHandle<edm4eic::MCRecoParticleAssociationCollection> m_inputParticleAssociation{
    "inputParticleAssociations",
    Gaudi::DataHandle::Reader,
    this};
  mutable DataHandle<edm4eic::InclusiveKinematicsCollection> m_outputInclusiveKinematicsCollection{
    "outputInclusiveKinematics",
    Gaudi::DataHandle::Writer,
    this};

  Gaudi::Property<double> m_crossingAngle{this, "crossingAngle", -0.025 * Gaudi::Units::radian};

  SmartIF<IParticleSvc> m_pidSvc;
  double m_proton{0}, m_neutron{0}, m_electron{0};

public:
  InclusiveKinematicsJB(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {
    declareProperty("inputMCParticles", m_inputMCParticleCollection, "MCParticles");
    declareProperty("inputReconstructedParticles", m_inputParticleCollection, "ReconstructedParticles");
    declareProperty("inputParticleAssociations", m_inputParticleAssociation, "MCRecoParticleAssociation");
    declareProperty("outputInclusiveKinematics", m_outputInclusiveKinematicsCollection, "InclusiveKinematicsJB");
  }

  StatusCode initialize() override {
    if (Gaudi::Algorithm::initialize().isFailure())
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

  StatusCode execute(const EventContext&) const override {
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
    //  [&ef_coll](const auto& a){ return a.getSimID() == ef_coll[0].getObjectID().index; });
    auto ef_assoc = rcassoc.begin();
    for (; ef_assoc != rcassoc.end(); ++ef_assoc) {
      if (ef_assoc->getSimID() == (unsigned) ef_coll[0].getObjectID().index) {
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
    const auto ef_rc_id{ef_rc.getObjectID().index};

    // Loop over reconstructed particles to get all outgoing particles other than the scattered electron
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

    // Get boost to colinear frame
    auto boost = Jug::Base::Boost::determine_boost(ei, pi);

    for (const auto& p: rcparts) {
      // Get the scattered electron index and angle
      if (p.getObjectID().index == ef_rc_id) {

      // Sum over all particles other than scattered electron
      } else {
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
    auto ptsum = sqrt(pxsum*pxsum + pysum*pysum);

    // Sigma zero or negative
    if (sigma_h <= 0) {
      if (msgLevel(MSG::DEBUG)) {
        debug() << "Sigma zero or negative" << endmsg;
      }
      return StatusCode::SUCCESS;
    }

    // Calculate kinematic variables
    const auto y_jb = sigma_h / (2.*ei.energy());
    const auto Q2_jb = ptsum*ptsum / (1. - y_jb);
    const auto x_jb = Q2_jb / (4.*ei.energy()*pi.energy()*y_jb);
    const auto nu_jb = Q2_jb / (2.*m_proton*x_jb);
    const auto W_jb = sqrt(m_proton*m_proton + 2*m_proton*nu_jb - Q2_jb);
    auto kin = out_kinematics.create(x_jb, Q2_jb, W_jb, y_jb, nu_jb);
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
DECLARE_COMPONENT(InclusiveKinematicsJB)

} // namespace Jug::Reco
