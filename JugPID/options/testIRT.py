from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

from Configurables import PodioInput
from Configurables import Jug__PID__TestIRT as TestIRT

geo_service  = GeoSvc("GeoSvc", detectors=["../athena/athena.xml"])
podioevent = EICDataSvc("EventDataSvc", inputs=["../dd4/erich/out/sim_run.root"], OutputLevel=DEBUG) # FIXME
podioinput = PodioInput("PodioReader", collections=["mcparticles", "ERICHHits"], OutputLevel=DEBUG)

irtrec = TestIRT(
        "irt_rec",
        inputHitCollection="ERICHHits"
        )

# REFERENCE
#emcaldigi = EcalTungstenSamplingDigi("ecal_digi",
                                     #inputHitCollection="EcalBarrelHits",
                                     #outputHitCollection="DigiEcalBarrelHits",
                                     #inputEnergyUnit=units.GeV,
                                     #inputTimeUnit=units.ns,
                                     #energyResolutions=[0., 0.02, 0.],
                                     #dynamicRangeADC=700*units.keV,
                                     #pedestalSigma=40,
                                     #OutputLevel=DEBUG)

out = PodioOutput("out", filename="reco_testIRT.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
        TopAlg = [podioinput, irtrec, out],
        EvtSel = 'NONE',
        EvtMax   = 5,
        ExtSvc = [podioevent],
        OutputLevel=DEBUG
        )
