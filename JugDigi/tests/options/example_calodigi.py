from Gaudi.Configuration import *

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

geo_service  = GeoSvc("GeoSvc")#detectors=["topside/vertex_tracker.xml"])
podioevent   = EICDataSvc("EventDataSvc", inputs=["derp.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Digi__ExampleCaloDigi as ExampleCaloDigi
#from Configurables import Jug__Digi__ExampleCaloDigiFunc as ExampleCaloDigiFunc
from Configurables import Jug__Digi__UFSDTrackerDigi as UFSDTrackerDigi
from Configurables import Jug__Reco__TrackerHitReconstruction as TrackerHitReconstruction
podioinput = PodioInput("PodioReader", collections=["MCParticles","LAEC_PrShHits","LAEC_ShHits","FAEC_PrShHits","FAEC_ShHits","GEMTrackerHits"], OutputLevel=DEBUG)
caldigi = ExampleCaloDigi(inputHitCollection="FAEC_ShHits",outputHitCollection="RawFAECShowerHits")
ufsd_digi = UFSDTrackerDigi(inputHitCollection="GEMTrackerHits",outputHitCollection="GEMRawHits")
#caldigifunc = ExampleCaloDigiFunc(InputData="FAEC_ShHits",OutputData="DERP")
trackerhit = TrackerHitReconstruction(inputHitCollection="GEMRawHits",outputHitCollection="GEMTrackHits")

types = []

# this printout is useful to check that the type information is passed to python correctly
print("---------------------------------------\n")
print("---\n# List of input and output types by class")
for configurable in sorted([
        PodioInput, EICDataSvc, PodioOutput,
        ExampleCaloDigi,ExampleCaloDigi, UFSDTrackerDigi ],
                           key=lambda c: c.getType()):
    print("\"{}\":".format(configurable.getType()))
    props = configurable.getDefaultProperties()
    for propname, prop in sorted(props.items()):
        print(" prop name: {}".format(propname))
        if isinstance(prop, DataObjectHandleBase):
            types.append(prop.type())
            print("  {}: \"{}\"".format(propname, prop.type()))
print("---")


out = PodioOutput("out", filename="test.root")
out.outputCommands = ["keep *"]


ApplicationMgr(
    TopAlg = [podioinput, caldigi,ufsd_digi,trackerhit, out
              ],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

