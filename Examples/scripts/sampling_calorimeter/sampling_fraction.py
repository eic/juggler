import numpy as np
import pandas as pd
import ROOT
from ROOT import gROOT, gInterpreter
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
gInterpreter.AddIncludePath("/home/cpeng/apps/eic/include");
gInterpreter.AddIncludePath("/home/cpeng/apps/include");
gInterpreter.Declare('''
    #include "eicd/ClusterData.h"
    #include "eicd/CalorimeterHitCollection.h"
    double sum_energy(ROOT::VecOps::RVec<eic::CalorimeterHitData> &tcol) {
        double energy = 0;
        for (auto &t : tcol) {
            energy += t.energy;
        }
        return energy;
    }

    double max_energy(ROOT::VecOps::RVec<eic::ClusterData> &tcol) {
        double energy = 0;
        for (auto &t : tcol) {
            if (t.energy > energy) energy = t.energy;
        }
        return energy;
    }
        ''')

rdf = rdf.Define('fraction', 'max_energy(EcalBarrelClusters)/5000')\
         .Define('r', 'sqrt(EcalBarrelClusters.position.x*EcalBarrelClusters.position.x + EcalBarrelClusters.position.y*EcalBarrelClusters.position.y)')\
         .Define('z', 'EcalBarrelClusters.position.z')\
         .Define('angle', 'acos(z/sqrt(r*r + z*z))/M_PI*180.')\
         .Define('etot', 'sum_energy(RecoEcalBarrelHits)')\
         .Define('fraction2', 'etot/5000')\
         .Define('mc_energy', 'sqrt(mcparticles.pex*mcparticles.pex + mcparticles.pey*mcparticles.pey + mcparticles.pez*mcparticles.pez)')


hist = rdf.Histo1D(ROOT.RDF.TH1DModel('fraction', 'Sampling Fraction;Fraction;Counts', 200, 0.001, 0.1), 'fraction2')
c1 = ROOT.TCanvas('c1', '', 2560, 1440)
hist.Fit('gaus')
hist.Draw()
c1.SaveAs('sampling_fraction.png')

hist = rdf.Histo1D(ROOT.RDF.TH1DModel('energy', ';Energy (MeV);Counts', 100, 0, 100), 'etot')
c1 = ROOT.TCanvas('c1', '', 2560, 1440)
hist.Draw()
c1.SaveAs('sampling_energy.png')

c1 = ROOT.TCanvas('c1', '', 2560, 1440)
hist = rdf.Histo1D(ROOT.RDF.TH1DModel('angle', ';Angle (deg);Counts', 180, 0, 180), 'angle')
hist.Draw()
c1.SaveAs('sampling_angle.png')

c1 = ROOT.TCanvas('c1', '', 2560, 1440)
hist = rdf.Histo2D(ROOT.RDF.TH2DModel('2d', ';Fraction;Angle (deg)', 40, 0.01, 0.1, 60, 60, 120), 'fraction', 'angle')
hist.Draw('colz')
c1.SaveAs('sampling_2d.png')

exit(1)

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
