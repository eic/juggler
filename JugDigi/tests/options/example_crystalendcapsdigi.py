from Gaudi.Configuration import *

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

geo_service  = GeoSvc("GeoSvc")
podioevent   = EICDataSvc("EventDataSvc", inputs=["output_emcal_electrons.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Digi__CrystalEndcapsDigi as CrystalEndcapsDigi

podioinput = PodioInput("PodioReader", collections=["MCParticles","EcalHits"], OutputLevel=DEBUG)
emcaldigi = CrystalEndcapsDigi(inputHitCollection="EcalHits",outputHitCollection="RawDigiEcalHits")

out = PodioOutput("out", filename="digi_emcal_electrons.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, emcaldigi, out],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

