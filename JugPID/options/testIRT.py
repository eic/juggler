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
from Configurables import Jug__Reco__TestIRTAlgorithm as TestIRTAlgorithm

# S13660-3050AE-08 SiPM quantum efficiency [(wavelength [nm], q.e.)]
# Note: this is consistent with S13361-3050AE-08 (for eRICH)
# TODO: is this where we want these parameters?
qe_data = [
    (325, 0.04),
    (340, 0.10),
    (350, 0.20),
    (370, 0.30),
    (400, 0.35),
    (450, 0.40),
    (500, 0.38),
    (550, 0.35),
    (600, 0.27),
    (650, 0.20),
    (700, 0.15),
    (750, 0.12),
    (800, 0.08),
    (850, 0.06),
    (900, 0.04)
]

podioinput = PodioInput(
        "PodioReader",
        collections=["mcparticles", "ERICHHits"],
        OutputLevel=DEBUG
        )

pmtdigi = PhotoMultiplierDigi(
        inputHitCollection="ERICHHits",
        outputHitCollection="DigiERICHHits",
        quantumEfficiency=[ ((1239.84/a)*units.eV, b) for a, b in qe_data ]
        )

pmtreco = PhotoMultiplierReco(
        inputHitCollection="DigiERICHHits",
        outputHitCollection="RecoERICHHits"
        )

irtrec = TestIRTAlgorithm(
        inputHitCollection="RecoERICHHits",
        outputClusterCollection="ERICHClusters"
        )

out = PodioOutput("out", filename=output_rec_file)
out.outputCommands = ["keep *"]

ApplicationMgr(
        TopAlg = [podioinput, pmtdigi, pmtreco, irtrec, out],
        EvtSel = 'NONE',
        EvtMax = n_events,
        ExtSvc = [podioevent],
        OutputLevel = DEBUG,
        PluginDebugLevel = 2
        )

