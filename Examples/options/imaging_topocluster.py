import os
import ROOT
from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc

from Configurables import PodioInput
from Configurables import Jug__Base__InputCopier_dd4pod__Geant4ParticleCollection_dd4pod__Geant4ParticleCollection_ as MCCopier
from Configurables import Jug__Base__InputCopier_dd4pod__CalorimeterHitCollection_dd4pod__CalorimeterHitCollection_ as CalCopier
from Configurables import Jug__Digi__CalorimeterHitDigi as CalorimeterHitDigi
from Configurables import Jug__Reco__ImagingPixelReco as ImagingPixelReco
from Configurables import Jug__Reco__ImagingTopoCluster as ImagingTopoCluster
from Configurables import Jug__Reco__ImagingClusterReco as ImagingClusterReco


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

print(kwargs)

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

imcaldigi = CalorimeterHitDigi("imcal_digi",
                               inputHitCollection="EcalBarrelHits",
                               outputHitCollection="DigiEcalBarrelHits",
                               inputEnergyUnit=units.GeV,
                               inputTimeUnit=units.ns,
                               energyResolutions=[0., 0.02, 0.],
                               dynamicRangeADC=3*units.MeV,
                               pedestalSigma=40,
                               OutputLevel=DEBUG)
imcalreco = ImagingPixelReco("imcal_reco",
                             inputHitCollection="DigiEcalBarrelHits",
                             outputHitCollection="RecoEcalBarrelHits",
                             dynamicRangeADC=3.*units.MeV,
                             pedestalSigma=40,
                             readoutClass="EcalBarrelHits",
                             layerField="layer",
                             sectorField="module")
imcalcluster = ImagingTopoCluster(inputHitCollection="RecoEcalBarrelHits",
                                  outputClusterCollection="EcalBarrelClusters",
                                  localRanges=[2.*units.mm, 2*units.mm],
                                  adjLayerRanges=[10*units.mrad, 10*units.mrad],
                                  adjLayerDiff=2,
                                  adjSectorDist=3.*units.cm)
clusterreco = ImagingClusterReco(inputClusterCollection="EcalBarrelClusters",
                                 outputLayerCollection="EcalBarrelClustersLayers",
                                 samplingFraction=sf,
                                 OutputLevel=DEBUG)

out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg=[podioinput, copier, calcopier, imcaldigi, imcalreco, imcalcluster, clusterreco, out],
    EvtSel='NONE',
    EvtMax=kwargs['nev'],
    ExtSvc=[podioevent],
    OutputLevel=DEBUG
)

