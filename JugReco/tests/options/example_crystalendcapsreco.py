from Gaudi.Configuration import *

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc
from GaudiKernel import SystemOfUnits as units

geo_service  = GeoSvc("GeoSvc", detectors=["../NPDet/src/GenericDetectors/calorimeters/compact/Crystal_example.xml"])
podioevent   = EICDataSvc("EventDataSvc", inputs=["output_emcal_electrons_npsim.root"], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Digi__CrystalEndcapsDigi as CrystalEndcapsDigi
from Configurables import Jug__Reco__CrystalEndcapsReco as CrystalEndcapsReco
from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster
from Configurables import Jug__Reco__ClusterRecoCoG as RecoCoG

podioinput = PodioInput("PodioReader", collections=["mcparticles","EcalHits"], OutputLevel=DEBUG)
emcaldigi = CrystalEndcapsDigi(inputHitCollection="EcalHits", outputHitCollection="RawDigiEcalHits")
emcalreco = CrystalEndcapsReco(inputHitCollection="RawDigiEcalHits", outputHitCollection="RecoEcalHits",
                               minModuleEdep=1.0*units.MeV)
emcalcluster = IslandCluster(inputHitCollection="RecoEcalHits", outputClusterCollection="EcalClusters",
                             minClusterCenterEdep=30*units.MeV, groupRange=2.0)
clusterreco = RecoCoG(clusterCollection="EcalClusters", logWeightThres=4.2)

out = PodioOutput("out", filename="reco_emcal_electrons_npsim.root")
out.outputCommands = ["keep EcalClusters"]

ApplicationMgr(
    TopAlg = [podioinput, emcaldigi, emcalreco, emcalcluster, clusterreco, out],
    EvtSel = 'NONE',
    EvtMax   = 100,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

