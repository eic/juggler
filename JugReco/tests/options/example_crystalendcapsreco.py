from Gaudi.Configuration import *

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

geo_service  = GeoSvc("GeoSvc", detectors=["../NPDet/src/GenericDetectors/calorimeters/compact/CrystalEndcapECAL_example.xml"])
podioevent   = EICDataSvc("EventDataSvc", inputs=["digi_emcal_electrons_npsim.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster

podioinput = PodioInput("PodioReader", collections=["RawDigiEcalHits"], OutputLevel=DEBUG)
emcalcluster = IslandCluster(inputHitCollection="RawDigiEcalHits", outputClusterCollection="EcalClusters")

out = PodioOutput("out", filename="reco_emcal_electrons_npsim.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, emcalcluster, out],
    EvtSel = 'NONE',
    EvtMax   = 100,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

