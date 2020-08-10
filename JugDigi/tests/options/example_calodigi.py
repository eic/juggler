from Gaudi.Configuration import *

from Configurables import ApplicationMgr, EICDataSvc, PodioOutput

podioevent   = EICDataSvc("EventDataSvc", inputs=["derp.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Digi__ExampleCaloDigi as ExampleCaloDigi
podioinput = PodioInput("PodioReader", collections=["mcparticles","FAEC_ShHits"], OutputLevel=DEBUG)
caldigi = ExampleCaloDigi(InputData="FAEC_ShHits",OutputData="RawFAECShowerHits")

out = PodioOutput("out", filename="test.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, caldigi, out
              ],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

