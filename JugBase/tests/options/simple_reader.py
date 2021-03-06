from Gaudi.Configuration import *
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput

podioevent   = EICDataSvc("EventDataSvc", inputs=["derp.root"], OutputLevel=DEBUG)

from Configurables import PodioInput, ReadTestConsumer
podioinput = PodioInput("PodioReader", collections=["mcparticles"], OutputLevel=DEBUG)
checker = ReadTestConsumer()

out = PodioOutput("out", filename="test.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, checker,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

