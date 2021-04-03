from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

from Configurables import PodioInput
from Configurables import Jug__Digi__EcalTungstenSamplingDigi as EcalTungstenSamplingDigi
from Configurables import Jug__Reco__EcalTungstenSamplingReco as EcalTungstenSamplingReco
from Configurables import Jug__Reco__SamplingECalHitsMerger as SamplingECalHitsMerger
from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster
from Configurables import Jug__Reco__ClusterRecoCoG as RecoCoG

geo_service  = GeoSvc("GeoSvc", detectors=["../topside/test.xml"])
podioevent = EICDataSvc("EventDataSvc", inputs=["barrel_electrons.root"], OutputLevel=DEBUG)

podioinput = PodioInput("PodioReader", collections=["mcparticles", "EcalBarrelHits"], OutputLevel=DEBUG)
emcaldigi = EcalTungstenSamplingDigi("ecal_digi",
                                     inputHitCollection="EcalBarrelHits",
                                     outputHitCollection="DigiEcalBarrelHits",
                                     inputEnergyUnit=units.GeV,
                                     inputTimeUnit=units.ns,
                                     OutputLevel=DEBUG)
emcalreco = EcalTungstenSamplingReco("ecal_reco",
                                     inputHitCollection="DigiEcalBarrelHits",
                                     outputHitCollection="RecoEcalBarrelHits",
                                     OutputLevel=DEBUG)
# readout id definition for barrel ecal
# <id>system:8,barrel:3,module:4,layer:6,slice:5,x:32:-16,y:-16</id>
# xy_merger sum layers/slices, masking (8+3+4, 8+3+4+5+6-1)
xymerger = SamplingECalHitsMerger("ecal_xy_merger",
                                   cellIDMaskRanges=[(15, 25)],
                                   inputHitCollection="RecoEcalBarrelHits",
                                   outputHitCollection="RecoEcalBarrelHitsXY")
# xy_merger sum modules, masking (8+3+4+5+6, 8+3+4+5+6+32-1)
zmerger = SamplingECalHitsMerger("ecal_z_merger",
                                  cellIDMaskRanges=[(26, 57)],
                                  inputHitCollection="RecoEcalBarrelHits",
                                  outputHitCollection="RecoEcalBarrelHitsZ")
emcalcluster = IslandCluster(inputHitCollection="RecoEcalBarrelHitsXY",
                             outputClusterCollection="EcalBarrelClusters",
                             minClusterCenterEdep=5.0*units.MeV,
                             splitCluster=False,
                             groupRange=5.0)
clusterreco = RecoCoG(clusterCollection="EcalBarrelClusters", logWeightBase=6.2, OutputLevel=DEBUG)


out = PodioOutput("out", filename="barrel_cluster.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, emcaldigi, emcalreco, xymerger, zmerger, emcalcluster, clusterreco, out],
    EvtSel = 'NONE',
    EvtMax   = 1000,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

