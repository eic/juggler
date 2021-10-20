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

geo_service  = GeoSvc("GeoSvc", detectors=["{}.xml".format(detector_name)])
podioevent   = EICDataSvc("EventDataSvc", inputs=[input_sim_file], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__PID__TestIRT as TestIRT

podioinput = PodioInput(
        "PodioReader",
        collections=["mcparticles", "ERICHHits"],
        OutputLevel=DEBUG
        )

irtrec = TestIRT(
        "irt_rec"#,
        #inputHitCollection="ERICHHits"
        )

out = PodioOutput("out", filename="reco_testIRT.root")
out.outputCommands = ["keep *"]
ApplicationMgr(
        TopAlg = [podioinput, irtrec, out],
        EvtSel = 'NONE',
        EvtMax   = n_events,
        ExtSvc = [podioevent],
        OutputLevel=DEBUG
        )

# REFERENCE: from reconstruction benchmark options file `rich_reco.py` #############
#
# from Configurables import Jug__Digi__PhotoMultiplierDigi as PhotoMultiplierDigi
# from Configurables import Jug__Reco__PhotoMultiplierReco as PhotoMultiplierReco
# from Configurables import Jug__Reco__PhotoRingClusters as PhotoRingClusters
#
# qe_data = [(1.0, 0.25), (7.5, 0.25),]
#
# pmtdigi = PhotoMultiplierDigi(inputHitCollection="ERICHHits", outputHitCollection="DigiForwardRICHHits",
#                               quantumEfficiency=[(a*units.eV, b) for a, b in qe_data])
# pmtreco = PhotoMultiplierReco(inputHitCollection="DigiForwardRICHHits", outputHitCollection="RecoForwardRICHHits")
# richcluster = PhotoRingClusters(inputHitCollection="RecoForwardRICHHits", #inputTrackCollection="mcparticles",
#                                 outputClusterCollection="RICHClusters")
# 
# out = PodioOutput("out", filename=output_rec_file)
# out.outputCommands = ["keep *"]
# 
# ApplicationMgr(
#     TopAlg = [podioinput, pmtdigi, pmtreco, richcluster, out],
#     EvtSel = 'NONE',
#     EvtMax = n_events,
#     ExtSvc = [podioevent],
#     OutputLevel=DEBUG
#  )

