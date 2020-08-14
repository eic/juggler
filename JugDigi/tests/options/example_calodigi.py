from Gaudi.Configuration import *

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput


podioevent   = EICDataSvc("EventDataSvc", inputs=["derp.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Digi__ExampleCaloDigi as ExampleCaloDigi
from Configurables import Jug__Digi__ExampleCaloDigiFunc as ExampleCaloDigiFunc
from Configurables import Jug__Digi__UFSDTrackerDigi as UFSDTrackerDigi
podioinput = PodioInput("PodioReader", collections=["mcparticles","LAEC_PrShHits","LAEC_ShHits","FAEC_PrShHits","FAEC_ShHits","GEMTrackerHits"], OutputLevel=DEBUG)
#caldigi = ExampleCaloDigi(inputHitCollection="FAEC_ShHits",outputHitCollection="RawFAECShowerHits")
ufsd_digi = UFSDTrackerDigi(inputHitCollection="GEMTrackerHits",outputHitCollection="GEMRawHits")
caldigifunc = ExampleCaloDigiFunc(InputData="FAEC_ShHits",OutputData="DERP")

types = []

# this printout is useful to check that the type information is passed to python correctly
print("---------------------------------------\n")
print("---\n# List of input and output types by class")
for configurable in sorted([
        PodioInput, EICDataSvc, PodioOutput,
        ExampleCaloDigiFunc,ExampleCaloDigi, UFSDTrackerDigi ],
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
    TopAlg = [podioinput, ufsd_digi,caldigifunc, out
              ],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

