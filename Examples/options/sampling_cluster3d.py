import os
from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

from Configurables import PodioInput
from Configurables import Jug__Base__InputCopier_dd4pod__Geant4ParticleCollection_dd4pod__Geant4ParticleCollection_ as MCCopier
from Configurables import Jug__Digi__EcalTungstenSamplingDigi as EcalTungstenSamplingDigi
from Configurables import Jug__Reco__EcalTungstenSamplingReco as EcalTungstenSamplingReco
from Configurables import Jug__Reco__TopologicalCellCluster as TopologicalCellCluster
from Configurables import Jug__Reco__ImagingClusterReco as ImagingReco

# get sampling fraction from system environment variable, 1.0 by default
sf = float(os.environ.get('CB_EMCAL_SAMP_FRAC', '1.0'))

geo_service  = GeoSvc("GeoSvc", detectors=["../topside/test.xml"])
podioevent = EICDataSvc("EventDataSvc", inputs=["../topside/test.root"], OutputLevel=DEBUG)
out = PodioOutput("out", filename="barrel_cluster.root")

podioinput = PodioInput("PodioReader", collections=["mcparticles", "EcalBarrelHits"], OutputLevel=DEBUG)

copier = MCCopier("MCCopier",
         inputCollection="mcparticles",
         outputCollection="mcparticles2",
         OutputLevel=DEBUG)
emcaldigi = EcalTungstenSamplingDigi("ecal_digi",
                                     inputHitCollection="EcalBarrelHits",
                                     outputHitCollection="DigiEcalBarrelHits",
                                     inputEnergyUnit=units.GeV,
                                     inputTimeUnit=units.ns,
                                     energyResolutions=[0., 0.02, 0.],
                                     dynamicRangeADC=3*units.MeV,
                                     pedestalSigma=40,
                                     OutputLevel=DEBUG)
emcalreco = EcalTungstenSamplingReco("ecal_reco",
                                     inputHitCollection="DigiEcalBarrelHits",
                                     outputHitCollection="RecoEcalBarrelHits",
                                     dynamicRangeADC=3.*units.MeV,
                                     pedestalSigma=40,
                                     OutputLevel=DEBUG)
emcalcluster = TopologicalCellCluster(inputHitCollection="RecoEcalBarrelHits",
                                      outputClusterCollection="EcalBarrelClusters",
                                      minClusterCenterEdep=0.5*units.MeV,
                                      groupRanges=[5.*units.cm, 5*units.cm, 5.*units.cm])
clusterreco = ImagingReco(inputClusterCollection="EcalBarrelClusters",
                          outputClusterCollection="EcalBarrelClustersReco",
                          outputLayerCollection="EcalBarrelClustersLayers",
                          samplingFraction=sf,
                          layerIDMaskRange=[15, 24],
                          OutputLevel=DEBUG)

out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, copier, emcaldigi, emcalreco, emcalcluster, clusterreco, out],
    EvtSel = 'NONE',
    EvtMax   = 10000,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )

