import os
from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

from Configurables import PodioInput
from Configurables import Jug__Digi__EcalTungstenSamplingDigi as EcalTungstenSamplingDigi
from Configurables import Jug__Reco__EcalTungstenSamplingReco as EcalTungstenSamplingReco
from Configurables import Jug__Reco__CalorimeterHitsMerger as CalorimeterHitsMerger
from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster
from Configurables import Jug__Reco__ClusterRecoCoG as RecoCoG

# get sampling fraction from system environment variable, 1.0 by default
sf = float(os.environ.get('CB_EMCAL_SAMP_FRAC', '1.0'))

geo_service  = GeoSvc("GeoSvc", detectors=["../athena/athena.xml"])
podioevent = EICDataSvc("EventDataSvc", inputs=["../reconstruction_benchmarks/sim_barrel_clusters.root"], OutputLevel=DEBUG)

podioinput = PodioInput("PodioReader", collections=["mcparticles", "EcalBarrelHits"], OutputLevel=DEBUG)
emcaldigi = EcalTungstenSamplingDigi("ecal_digi",
                                     inputHitCollection="EcalBarrelHits",
                                     outputHitCollection="DigiEcalBarrelHits",
                                     inputEnergyUnit=units.GeV,
                                     inputTimeUnit=units.ns,
                                     energyResolutions=[0., 0.02, 0.],
                                     dynamicRangeADC=700*units.keV,
                                     pedestalSigma=40,
                                     OutputLevel=DEBUG)
emcalreco = EcalTungstenSamplingReco("ecal_reco",
                                     inputHitCollection="DigiEcalBarrelHits",
                                     outputHitCollection="RecoEcalBarrelHits",
                                     dynamicRangeADC=700*units.keV,
                                     pedestalSigma=40,
                                     OutputLevel=DEBUG)
# readout id definition for barrel ecal
# <id>system:8,barrel:3,module:4,layer:10,slice:5,x:32:-16,y:-16</id>
# sum layers/slices, masking (8+3+4, 8+3+4+5+10-1)
xymerger = CalorimeterHitsMerger("ecal_xy_merger",
                                 fields=["layer", "slice"],
                                 inputHitCollection="RecoEcalBarrelHits",
                                 outputHitCollection="RecoEcalBarrelHitsXY")
# sum modules, masking (8+3+4+5+10, 8+3+4+5+10+32-1)
zmerger = CalorimeterHitsMerger("ecal_z_merger",
                                fields=["x", "y"],
                                inputHitCollection="RecoEcalBarrelHits",
                                outputHitCollection="RecoEcalBarrelHitsZ")
emcalcluster = IslandCluster(inputHitCollection="RecoEcalBarrelHitsXY",
                             outputClusterCollection="EcalBarrelClusters",
                             minClusterCenterEdep=0.5*units.MeV,
                             splitCluster=False,
                             groupRanges=[5.*units.cm, 5*units.cm, 5.*units.cm])
clusterreco = RecoCoG(clusterCollection="EcalBarrelClusters", logWeightBase=6.2, samplingFraction=sf, OutputLevel=DEBUG)


out = PodioOutput("out", filename="barrel_cluster.root")
out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, emcaldigi, emcalreco, xymerger, zmerger, emcalcluster, clusterreco, out],
    EvtSel = 'NONE',
    EvtMax   = 5,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

