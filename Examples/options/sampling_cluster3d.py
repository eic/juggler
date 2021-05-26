import os
import ROOT
from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

from Configurables import PodioInput
from Configurables import Jug__Base__InputCopier_dd4pod__Geant4ParticleCollection_dd4pod__Geant4ParticleCollection_ as MCCopier
from Configurables import Jug__Base__InputCopier_dd4pod__CalorimeterHitCollection_dd4pod__CalorimeterHitCollection_ as CalCopier
from Configurables import Jug__Digi__EcalTungstenSamplingDigi as EcalTungstenSamplingDigi
from Configurables import Jug__Reco__EcalTungstenSamplingReco as EcalTungstenSamplingReco
from Configurables import Jug__Reco__TopologicalCellCluster as TopologicalCellCluster
from Configurables import Jug__Reco__ImagingClusterReco as ImagingReco


# input arguments through environment variables
kwargs = dict()
kwargs['sf'] = float(os.environ.get('CB_EMCAL_SAMP_FRAC', '1.0'))
kwargs['input'] = os.environ.get('CB_EMCAL_SIM_FILE', '../topside/barrel_pion0_5GeV.root')
kwargs['output'] = os.environ.get('CB_EMCAL_REC_FILE', 'barrel_pion0_5GeV_cluster.root')
kwargs['compact'] = os.environ.get('CB_EMCAL_COMPACT_PATH', '../topside/test.xml')
kwargs['nev'] = int(os.environ.get('CB_EMCAL_NUMEV', 10000))

if kwargs['nev'] < 1:
    f = ROOT.TFile(kwargs['input'])
    kwargs['nev'] = f.events.GetEntries()

# get sampling fraction from system environment variable, 1.0 by default
sf = float(os.environ.get('CB_EMCAL_SAMP_FRAC', '1.0'))

geo_service = GeoSvc("GeoSvc", detectors=kwargs['compact'].split(','))
podioevent = EICDataSvc("EventDataSvc", inputs=kwargs['input'].split(','), OutputLevel=DEBUG)
out = PodioOutput("out", filename=kwargs['output'])

podioinput = PodioInput("PodioReader", collections=["mcparticles", "EcalBarrelHits"], OutputLevel=DEBUG)

copier = MCCopier("MCCopier",
                  inputCollection="mcparticles",
                  outputCollection="mcparticles2",
                  OutputLevel=DEBUG)
calcopier = CalCopier("CalCopier",
                      inputCollection="EcalBarrelHits",
                      outputCollection="EcalBarrelHits2",
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
                                      minClusterCenterEdep=0.3*units.MeV,
                                      localRanges=[2.*units.mm, 2*units.mm],
                                      adjLayerRanges=[5*units.mrad, 5*units.mrad],
                                      adjLayerDiff=1,
                                      adjSectorDist=1*units.cm,
                                      layerField="layer",
                                      sectorField="module")
clusterreco = ImagingReco(inputClusterCollection="EcalBarrelClusters",
                          outputClusterCollection="EcalBarrelClustersReco",
                          outputLayerCollection="EcalBarrelClustersLayers",
                          samplingFraction=sf,
                          layerIDMaskRange=[15, 24],
                          OutputLevel=DEBUG)

out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg=[podioinput, copier, calcopier, emcaldigi, emcalreco, emcalcluster, clusterreco, out],
    EvtSel='NONE',
    EvtMax=kwargs['nev'],
    ExtSvc=[podioevent],
    OutputLevel=DEBUG
)

