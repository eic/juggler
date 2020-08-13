from Gaudi.Configuration import *

from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

podioevent   = EICDataSvc("EventDataSvc", inputs=["derp.root"], OutputLevel=DEBUG)
geo_service  = GeoSvc(detectors=["topside/vertex_tracker.xml"])

# reads HepMC text file and write the HepMC::GenEvent to the data service
from Configurables import PodioInput, ReadTestConsumer
podioinput = PodioInput("PodioReader", collections=["mcparticles","FAEC_ShHits"], OutputLevel=DEBUG)
checker = ReadTestConsumer()

out = PodioOutput("out", filename="test.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, checker,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent,geo_service],
    OutputLevel=DEBUG
 )

