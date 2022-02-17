from Gaudi.Configuration import *

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

from Configurables import PodioInput
from Configurables import Jug__Digi__CrystalEndcapsDigi as CrystalEndcapsDigi

#geo_service  = GeoSvc("GeoSvc")
podioevent   = EICDataSvc("EventDataSvc", inputs=["output_emcal_electrons_npsim.root"], OutputLevel=DEBUG)

podioinput = PodioInput("PodioReader", collections=["MCParticles","EcalHits"], OutputLevel=DEBUG)
emcaldigi = CrystalEndcapsDigi("ecal_digi",inputHitCollection="EcalHits",outputHitCollection="RawDigiEcalHits", OutputLevel=DEBUG)

out = PodioOutput("out", filename="digi_emcal_electrons_npsim.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, emcaldigi, out],
    EvtSel = 'NONE',
    EvtMax   = 100,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

