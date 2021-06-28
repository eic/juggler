// A specialized digitization for GlueX-like barrel ecal
//
// Author: Chao Peng, Maria Zurek (ANL)
// Date: 06/27/2021

#include <algorithm>
#include <cmath>

#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/RndmGenerators.h"

// FCCSW
#include "JugBase/DataHandle.h"
#include "JugBase/IGeoSvc.h"

// Event Model related classes
#include "dd4pod/CalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitData.h"

using namespace Gaudi::Units;

namespace Jug::Digi {

class ScFiBarrelCalDigi : public GaudiAlgorithm {
public:
    // additional smearing resolutions
    Gaudi::Property<std::vector<double>>        u_eRes{this, "energyResolutions", {}}; // a/sqrt(E/GeV) + b + c/(E/GeV)
    Gaudi::Property<double>                     m_tRes{this, "timingResolution", 0.0*ns};

    // @TODO, this may be a vector of parameters, need some reference
    Gaudi::Property<std::vector<double>>        u_atten{this, "lightAttenuation", {0.}};
    // geometry service to decode readout ids, merging fibers with the same light guide
    Gaudi::Property<std::string>                m_geoSvcName{this, "geoServiceName", "GeoSvc"};
    Gaudi::Property<std::string>                m_readout{this, "readoutClass", ""};
    Gaudi::Property<std::vector<std::string>>   u_mfields{this, "mergingFields", {}};
    SmartIF<IGeoSvc> m_geoSvc;

    // input units, should not be changed
    Gaudi::Property<double>                     m_eUnit{this, "inputEnergyUnit", GeV};
    Gaudi::Property<double>                     m_tUnit{this, "inputTimeUnit", ns};

    // digitization settings
    Gaudi::Property<int>                        m_capADC{this, "capacityADC", 8096};
    Gaudi::Property<double>                     m_dyRangeADC{this, "dynamicRangeADC", 100*MeV};
    Gaudi::Property<int>                        m_pedMeanADC{this, "pedestalMean", 400};
    Gaudi::Property<double>                     m_pedSigmaADC{this, "pedestalSigma", 3.2};
    Rndm::Numbers                               m_normDist;

    DataHandle<dd4pod::CalorimeterHitCollection>
        m_inputHitCollection{"inputHitCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<eic::RawCalorimeterHitCollection>
        m_outputHitCollection{"outputHitCollection", Gaudi::DataHandle::Writer, this};
    double res[3] = {0., 0., 0.};

      //  ill-formed: using GaudiAlgorithm::GaudiAlgorithm;
    ScFiBarrelCalDigi(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc)
    {
        declareProperty("inputHitCollection", m_inputHitCollection, "");
        declareProperty("outputHitCollection", m_outputHitCollection, "");
    }

    StatusCode initialize() override
    {
        if (GaudiAlgorithm::initialize().isFailure()) {
            return StatusCode::FAILURE;
        }
        // random number generator from service
        auto randSvc = svc<IRndmGenSvc>("RndmGenSvc", true);
        auto sc = m_normDist.initialize(randSvc, Rndm::Gauss(0.0, 1.0));
        if (!sc.isSuccess()) {
            return StatusCode::FAILURE;
        }
        // set energy resolution numbers
        for (size_t i = 0; i < u_eRes.size() && i < 3; ++i) {
            res[i] = u_eRes[i];
        }

        // geometry service
        m_geoSvc = service(m_geoSvcName);
        if (!m_geoSvc) {
            error() << "Unable to locate Geometry Service. "
                    << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
            return StatusCode::FAILURE;
        }

        if (m_readout.value().empty() || u_mfields.value().empty()) {
            warning() << "No readout class or ID fields provided, taking single fiber signals without any merge."
                      << endmsg;
            return StatusCode::SUCCESS;
        }

        // @TODO: get id mask
        return StatusCode::SUCCESS;
    }

    StatusCode execute() override
    {
        // @TODO: implement light attenuation and grid merging
        // input collections
        const auto simhits = m_inputHitCollection.get();
        // Create output collections
        auto rawhits = m_outputHitCollection.createAndPut();
        for (const auto& ahit : *simhits) {
            double eres = std::sqrt(std::pow(m_normDist()*res[0] / sqrt(ahit.energyDeposit()*m_eUnit/GeV), 2)
                                    + std::pow(m_normDist()*res[1], 2)
                                    + std::pow(m_normDist()*res[2] / (ahit.energyDeposit()*m_eUnit/GeV), 2));
            double ped = m_pedMeanADC + m_normDist()*m_pedSigmaADC;
            long long adc = std::llround(ped + ahit.energyDeposit()*(1. + eres) * m_eUnit/m_dyRangeADC*m_capADC);
            eic::RawCalorimeterHit rawhit(
                (long long) ahit.cellID(),
                (adc > m_capADC ? m_capADC.value() : adc),
                (double) ahit.truth().time*m_tUnit/ns + m_normDist()*m_tRes/ns
                );
            rawhits->push_back(rawhit);
        }
        return StatusCode::SUCCESS;
    }
};

DECLARE_COMPONENT(ScFiBarrelCalDigi)

} // namespace Jug::Digi
