from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

detector_name = "athena"
if "JUGGLER_DETECTOR" in os.environ:
  detector_name = str(os.environ["JUGGLER_DETECTOR"])
if "DETECTOR_PATH" in os.environ:
  detector_name = str(os.environ["DETECTOR_PATH"])+"/"+detector_name

# TODO: check for arguments, and catch errors on unset env vars
input_sim_file  = str(os.environ["JUGGLER_SIM_FILE"])
output_rec_file = str(os.environ["JUGGLER_REC_FILE"])
n_events = str(os.environ["JUGGLER_N_EVENTS"])

geo_service = GeoSvc("GeoSvc", detectors=["{}.xml".format(detector_name)])
podioevent = EICDataSvc("EventDataSvc", inputs=[input_sim_file], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Digi__PhotoMultiplierDigi as PhotoMultiplierDigi
from Configurables import Jug__Reco__PhotoMultiplierReco as PhotoMultiplierReco
from Configurables import Jug__Reco__TestIRT as TestIRT

podioinput = PodioInput(
        "PodioReader",
        collections=["mcparticles", "ERICHHits"],
        OutputLevel=DEBUG
        )

irtrec = TestIRT(
        "irt_rec",
        inputHitCollection="ERICHHits",
        outputPidCollection="ERICHPid"
        )

out = PodioOutput("out", filename=output_rec_file)
out.outputCommands = ["keep *"]

ApplicationMgr(
        TopAlg = [podioinput, irtrec, out],
        EvtSel = 'NONE',
        EvtMax = n_events,
        ExtSvc = [podioevent],
        OutputLevel = DEBUG
        )

