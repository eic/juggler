from Gaudi.Configuration import *
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput

podioevent   = EICDataSvc("EventDataSvc", inputs=["derp.root"], OutputLevel=DEBUG)

from Configurables import Jug__Base__InputCopier_dd4pod__Geant4ParticleCollection_dd4pod__Geant4ParticleCollection_ as MCCopier
from Configurables import PodioInput, ReadTestConsumer
podioinput = PodioInput("PodioReader", collections=["mcparticles"], OutputLevel=DEBUG)
checker = ReadTestConsumer()
copier = MCCopier("copier", inputCollection="mcparticles", outputCollection="mcparticles2",OutputLevel=DEBUG)

out = PodioOutput("out", filename="test.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, checker,copier,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

