import numpy as np
import pandas as pd
import ROOT
from ROOT import gROOT
import argparse
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

parser = argparse.ArgumentParser(description='sampling calorimeter performance')
parser.add_argument('file', type=str, help='path to root file')
args = parser.parse_args()

gROOT.Macro('rootlogon.C')
# root dataframe
rdf = ROOT.RDataFrame('events', args.file)

rdf = rdf.Define('fraction', 'EcalBarrelClusters.energy/5000')\
         .Define('r', 'sqrt(EcalBarrelClusters.position.x*EcalBarrelClusters.position.x + EcalBarrelClusters.position.y*EcalBarrelClusters.position.y)')\
         .Define('z', 'EcalBarrelClusters.position.z')\
         .Define('angle', 'acos(z/sqrt(r*r + z*z))/M_PI*180.')

hist = rdf.Histo1D(ROOT.RDF.TH1DModel('energy', 'Sampling Fraction;Fraction;Counts', 280, 0.04, 0.6), 'fraction')
c1 = ROOT.TCanvas('c1', '', 2560, 1440)
hist.Fit('gaus')
hist.Draw()
c1.SaveAs('sampling_energy.png')

c1 = ROOT.TCanvas('c1', '', 2560, 1440)
hist = rdf.Histo1D(ROOT.RDF.TH1DModel('angle', ';Angle (deg);Counts', 180, 0, 180), 'angle')
hist.Draw()
c1.SaveAs('sampling_angle.png')

c1 = ROOT.TCanvas('c1', '', 2560, 1440)
hist = rdf.Histo2D(ROOT.RDF.TH2DModel('2d', ';Fraction;Angle (deg)', 40, 0.1, 0.5, 60, 60, 120), 'fraction', 'angle')
hist.Draw('colz')
c1.SaveAs('sampling_2d.png')

f = ROOT.TFile(args.file)
events = f.events

for iev in range(events.GetEntries()):
    events.GetEntry(iev)
    xedep = []
    for hit in events.RecoEcalBarrelHitsZ:
        xedep.append((np.sqrt(hit.position.x**2 + hit.position.y**2), hit.energy))
    df = pd.DataFrame(data=xedep, columns=['x', 'edep']).set_index('x', drop=True)
    df.sort_index(inplace=True)
    df.index -= 890
    print(df.cumsum())

fig, ax = plt.subplots()
ax.plot(df.index, df.cumsum()['edep'].values)
fig.savefig('test.png')
