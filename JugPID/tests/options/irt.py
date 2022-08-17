import os
from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc   
from Configurables import PodioInput
from Configurables import Jug__PID__IRTAlgorithm as IRTAlgorithm 

### SETTINGS ####################
compact_file = os.environ['DETECTOR_PATH'] + '/' + os.environ['DETECTOR'] +'_drich_only.xml'
sim_file     = os.environ['JUGGLER_SIM_FILE']
irt_auxfile  = os.environ['JUGGLER_IRT_AUXFILE']
out_file     = os.environ['JUGGLER_REC_FILE']
#################################

# Well, for now only need DRICH geometry for this example;
geo_service = GeoSvc("GeoSvc", detectors=[compact_file])

# Input file after 'npsim' pass;
podioevent = EICDataSvc("EventDataSvc", inputs=[sim_file])

# S13660-3050AE-08 SiPM quantum efficiency [(wavelength [nm], q.e.)]
# Note: this is consistent with S13361-3050AE-08 (for DRICH)
qe_data = [
    (325, 0.04),
    (340, 0.10),
    (350, 0.20),
    (370, 0.30),
    (400, 0.35),
    (450, 0.40),
    (500, 0.38),
    (550, 0.35),
    (600, 0.27),
    (650, 0.20),
    (700, 0.15),
    (750, 0.12),
    (800, 0.08),
    (850, 0.06),
    (900, 0.04)
]

radiators = [
    "Aerogel   zbins=5  smearing=gaussian 2mrad rindex=1.0190 attenuation[mm]=48.0",
    "GasVolume zbins=10 smearing=gaussian 5mrad rindex=1.00076"
]

podioinput = PodioInput(
        "PodioReader",
        # Input collections: MC truth tracks and DRICH raw hits (photons);
        collections=["MCParticles", "DRICHHits"],
        OutputLevel=DEBUG
        )

irtrec = IRTAlgorithm(
        # Input collections: MC truth tracks and DRICH raw hits (photons);
        inputMCParticles="MCParticles",

        # DRICH optics configuration produced by DRICH_geo.cpp code along with the dd4hep XML file;
        ConfigFile=irt_auxfile,
        Radiators=[ (r) for r in radiators ],

        # SiPM PDE; FIXME: units.eV coefficient gives extra x1000 (?);
        QEcurve=[ ((1239.84/a), b) for a, b in qe_data ],
        # Rebin the QE in that many equidistant bins internally;
        QEbins="100",
        # SiPM geometric fill factor and "safety factor" for the photon count estimates;
        #GeometricEfficiency="1.00",
        SafetyFactor="0.70",
        )

# Output ROOT file; keep the input collections as well, append DRICH PID tables;
out = PodioOutput(
        "out",
        filename=out_file,
        OutputLevel=DEBUG
        )

out.outputCommands = ["keep *"]

ApplicationMgr(
        TopAlg = [podioinput, irtrec, out],
        EvtSel = 'NONE',
        # Process that many events;
        EvtMax = 50000,
        ExtSvc = [podioevent],
        OutputLevel = DEBUG,
        PluginDebugLevel = 2
        )

