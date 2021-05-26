import os
import ROOT
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
from Configurables import Jug__Reco__ImagingCalMLDataSorter as MLDataSorter


# input arguments through environment variables
kwargs = dict()
kwargs['input'] = os.environ.get('CB_EMCAL_SIM_FILE', '../topside/ml_test_pin_100MeV.root')
kwargs['output'] = os.environ.get('CB_EMCAL_REC_FILE', 'ml_test_pin_100MeV.root')
kwargs['compact'] = os.environ.get('CB_EMCAL_COMPACT_PATH', '../topside/test.xml')
kwargs['nev'] = int(os.environ.get('CB_EMCAL_NUMEV', 100000))

if kwargs['nev'] < 1:
    f = ROOT.TFile(kwargs['input'])
    kwargs['nev'] = f.events.GetEntries()

print(kwargs)
# get sampling fraction from system environment variable, 1.0 by default
sf = float(os.environ.get('CB_EMCAL_SAMP_FRAC', '1.0'))

geo_service  = GeoSvc("GeoSvc", detectors=kwargs['compact'].split(','))
podioevent = EICDataSvc("EventDataSvc", inputs=kwargs['input'].split(','), OutputLevel=DEBUG)
out = PodioOutput("out", filename=kwargs['output'])

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
emcaldata = MLDataSorter("ecal_ml",
                         inputHitCollection="RecoEcalBarrelHits",
                         outputHitCollection="EcalBarrelHitsML",
                         layerIDMaskRange=[15, 24],
                         etaSize=0.001,
                         phiSize=0.001,
                         OutputLevel=DEBUG)

out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg=[podioinput, copier, emcaldigi, emcalreco, emcaldata, out],
    EvtSel='NONE',
    EvtMax=kwargs['nev'],
    ExtSvc=[podioevent],
    OutputLevel=DEBUG
)

