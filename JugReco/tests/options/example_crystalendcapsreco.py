from Gaudi.Configuration import *

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

geo_service  = GeoSvc("GeoSvc", detectors=["../NPDet/src/GenericDetectors/calorimeters/compact/Crystal_example.xml"])
podioevent   = EICDataSvc("EventDataSvc", inputs=["output_emcal_electrons_npsim.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster
from Configurables import Jug__Digi__CrystalEndcapsDigi as CrystalEndcapsDigi

podioinput = PodioInput("PodioReader", collections=["mcparticles","EcalHits"], OutputLevel=DEBUG)
emcaldigi = CrystalEndcapsDigi(inputHitCollection="EcalHits", outputHitCollection="RawDigiEcalHits")
emcalcluster = IslandCluster(inputHitCollection="RawDigiEcalHits", outputClusterCollection="EcalClusters")

out = PodioOutput("out", filename="reco_emcal_electrons_npsim.root")
out.outputCommands = ["keep EcalClusters"]

ApplicationMgr(
    TopAlg = [podioinput, emcaldigi, emcalcluster, out],
    EvtSel = 'NONE',
    EvtMax   = 100,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

