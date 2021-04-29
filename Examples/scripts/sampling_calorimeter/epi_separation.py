import os
import numpy as np
import pandas as pd
import ROOT
from ROOT import gROOT, gInterpreter
import argparse
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, MaxNLocator
import sys


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


# define Ecal shapes
rmin, thickness, length = 890, 20*(10. + 1.65), 860*2+500

parser = argparse.ArgumentParser(description='epi_separation')
parser.add_argument('efile', type=str, help='path to root file (electrons)')
parser.add_argument('pifile', type=str, help='path to root file (pions)')
parser.add_argument('--e-prof', type=str, default='e_prof_5gev.csv', help='path to electron profile')
parser.add_argument('--pi-prof', type=str, default='pi_prof_5gev.csv', help='path to electron profile')
parser.add_argument('--plot-dir', type=str, default='./plots', dest='outdir', help='output directory')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)


def get_data(path):
    f = ROOT.TFile(path)
    events = f.events
    # prepare data
    df1 = pd.DataFrame()
    for iev in np.arange(events.GetEntries()):
        events.GetEntry(iev)
        data = []
        for hit in events.RecoEcalBarrelHits:
            lid = (hit.cellID >> 15) & 1023
            data.append((lid, hit.position.x, hit.position.y, hit.position.z, hit.energy))
        df2 = pd.DataFrame(data=data, columns=['id', 'x', 'y', 'z', 'edep'])
        dft = df2.groupby('id')['edep'].sum().to_frame()
        dft.loc[:, 'event'] = iev
        df1 = pd.concat([df1, dft])
        # mc_p = (events.mcparticles2[2].psx, events.mcparticles2[2].psy, events.mcparticles2[2].psz)
    df1.reset_index(inplace=True)
    return df1


profe = pd.read_csv(args.e_prof).set_index('id', drop=True)
profpi = pd.read_csv(args.pi_prof).set_index('id', drop=True)
dfe = get_data(args.efile)
dfpi = get_data(args.pifile)

colors = ['royalblue', 'indianred']

energy = 5000.
fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
for prof, color, title in zip([profe, profpi], colors, ['electron', 'pion']):
    ax.plot(prof.index, prof['mean']/energy*100., color=color, lw=2)
    ax.fill_between(prof.index, (prof['mean'] - prof['std'])/energy*100., (prof['mean'] + prof['std'])/energy*100.,
                    color=color, alpha=0.3, label=title)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.grid(linestyle=':', which='both')
ax.tick_params(labelsize=24)
ax.set_xlabel('Layer', fontsize=26)
ax.set_ylabel('Sampling Fraction (%)', fontsize=26)
ax.legend(fontsize=26, loc='upper left')
fig.savefig(os.path.join(args.outdir, 'energy_profile.png'))


def calc_chi2(grp, idmin=5):
    grp2 = grp[grp['id'] >= idmin]
    mean, std = profe.loc[grp2['id'], ['mean', 'std']].values.T
    edeps = grp2['edep'].values
    return np.sqrt(np.sum((edeps - mean)**2/std**2)/float(len(edeps)))


fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
chi2 = dfe.groupby('event').apply(calc_chi2)
ax.hist(chi2, bins=np.linspace(0, 5, 200), weights=[1/float(len(chi2))]*chi2,
        ec='royalblue', color='royalblue', alpha=0.5, label='electrons')
chi2 = dfpi.groupby('event').apply(calc_chi2)
print(chi2[chi2 < 0.7])
ax.hist(chi2, bins=np.linspace(0, 5, 200), weights=[1/float(len(chi2))]*chi2,
        ec='indianred', color='indianred', alpha=0.5, label='pions')
ax.grid(linestyle=':', which='major')
ax.set_axisbelow(True)
ax.tick_params(labelsize=24)
ax.set_xlabel(r'$\chi^2$', fontsize=26)
ax.set_ylabel('Normalized Counts', fontsize=26)
# ax.set_yscale('log')
# ax.set_ylim(1e-3, 1.)
ax.legend(fontsize=26)
fig.savefig(os.path.join(args.outdir, 'edep_chi2.png'))

