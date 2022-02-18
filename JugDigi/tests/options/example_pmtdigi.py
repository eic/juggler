from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

geo_service  = GeoSvc("GeoSvc")
podioevent   = EICDataSvc("EventDataSvc", inputs=["rich_test.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Digi__PhotoMultiplierDigi as PhotoMultiplierDigi

qe_data = [(1.0, 0.25), (7.5, 0.25),]
podioinput = PodioInput("PodioReader", collections=["MCParticles", "ForwardRICHHits"], OutputLevel=DEBUG)
pmtdigi = PhotoMultiplierDigi(inputHitCollection="ForwardRICHHits", outputHitCollection="DigiForwardRICHHits",
                              quantumEfficiency=[(a*units.eV, b) for a, b in qe_data])

out = PodioOutput("out", filename="digi_rich_test.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, pmtdigi, out],
    EvtSel = 'NONE',
    EvtMax   = 100,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

