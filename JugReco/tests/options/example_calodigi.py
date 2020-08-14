from Gaudi.Configuration import *

from Configurables import ApplicationMgr, EICDataSvc, PodioOutput

podioevent   = EICDataSvc("EventDataSvc", inputs=["derp.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Digi__ExampleCaloDigi as ExampleCaloDigi
from Configurables import Jug__Digi__UFSDTrackerDigi as UFSDTrackerDigi
podioinput = PodioInput("PodioReader", collections=["mcparticles","LAEC_PrShHits","LAEC_ShHits","FAEC_PrShHits","FAEC_ShHits","GEMTrackerHits"], OutputLevel=DEBUG)
caldigi = ExampleCaloDigi(inputHitCollection="FAEC_ShHits",outputHitCollection="RawFAECShowerHits")
ufsd_digi = UFSDTrackerDigi(inputHitCollection="GEMTrackerHits",outputHitCollection="GEMRawHits")

out = PodioOutput("out", filename="test.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, caldigi,ufsd_digi, out
              ],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

