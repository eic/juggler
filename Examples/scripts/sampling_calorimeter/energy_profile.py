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

parser = argparse.ArgumentParser(description='sampling calorimeter energy profiling')
parser.add_argument('file', type=str, help='path to root file')
parser.add_argument('--plot-dir', type=str, default='./plots', dest='outdir', help='output directory')
parser.add_argument('--save-prof', type=str, default='', dest='save', help='path to save profile')
parser.add_argument('--color', type=str, default='royalblue', dest='color', help='colors for bar plots')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

gROOT.Macro('rootlogon.C')
f = ROOT.TFile(args.file)
events = f.events

# prepare data
dfe = pd.DataFrame()
for iev in np.arange(events.GetEntries()):
    events.GetEntry(iev)
    data = []
    for hit in events.RecoEcalBarrelHits:
        lid = (hit.cellID >> 15) & 1023
        data.append((lid, hit.position.x, hit.position.y, hit.position.z, hit.energy))
    df = pd.DataFrame(data=data, columns=['id', 'x', 'y', 'z', 'edep'])
    dft = df.groupby('id')['edep'].sum().to_frame()
    dft.loc[:, 'event'] = iev
    dfe = pd.concat([dfe, dft])
    # mc_p = (events.mcparticles2[2].psx, events.mcparticles2[2].psy, events.mcparticles2[2].psz)
dfe.reset_index(inplace=True)
# print(dfe)


def find_start_layer(grp, min_edep=0.5):
    ids, edeps = grp.sort_values('id')[['id', 'edep']].values.T
    edeps = np.cumsum(edeps)
    return min(ids[edeps > min_edep]) if sum(edeps > min_edep) > 0 else 20


slayer = dfe.groupby('event').apply(find_start_layer).astype(int)
dfe = dfe.merge(slayer.to_frame(name='slayer'), on='event')
# dfe.loc[:, 'cid'] = dfe['id'] - dfe['slayer']
# prof = dfe[dfe['cid'] > 0].groupby('cid')['edep'].describe()
prof = dfe.groupby('id')['edep'].describe()
# print(prof['mean'])
bpos, bprops = (0.5, 0.95), dict(boxstyle='round', facecolor='wheat', alpha=0.3)


fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
ax.hist(dfe.groupby('event')['slayer'].min(), weights=[1/float(events.GetEntries())]*events.GetEntries(),
        ec='black', bins=np.arange(0.5, 10.5, step=1.0))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.grid(linestyle=':', which='both')
ax.tick_params(labelsize=24)
ax.set_xlabel('Start Layer', fontsize=26)
ax.set_ylabel('Normalized Counts', fontsize=26)
ax.text(*bpos, 'Mininum Edep\n' + '{:.1f} MeV'.format(1.0),
        transform=ax.transAxes, fontsize=26, verticalalignment='top', bbox=bprops)
fig.savefig(os.path.join(args.outdir, 'edep_start.png'))


fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
ax.plot(prof.index, prof['mean'], color=args.color, lw=2)
# ax.fill_between(prof.index, prof['25%'], prof['75%'], color=args.color, alpha=0.3)
ax.fill_between(prof.index, prof['mean'] - prof['std'], prof['mean'] + prof['std'], color=args.color, alpha=0.3)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.grid(linestyle=':', which='both')
ax.tick_params(labelsize=24)
ax.set_xlabel('Layer', fontsize=26)
ax.set_ylabel('Energy Deposit (MeV)', fontsize=26)
fig.savefig(os.path.join(args.outdir, 'edep.png'))

layers = np.asarray([
    [1, 5, 8,],
    [10, 15, 20],
])

fig, ax = plt.subplots(*layers.shape, figsize=(16, 9), dpi=160, sharex='col', sharey='all',
                       gridspec_kw=dict(hspace=0.05, wspace=0.05))

for ax, layer in zip(ax.flat, layers.flatten()):
    data = dfe[dfe['id'] == layer]
    ax.hist(data['edep'], weights=[1/float(len(data))]*len(data), bins=np.linspace(0, 20, 40),
            ec='black', color=args.color)
    ax.tick_params(labelsize=24)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.grid(linestyle=':', which='both')
    # ax.set_xlabel('Energy Deposit (MeV)', fontsize=26)
    # ax.set_ylabel('Normalized Counts', fontsize=26)
    mean, std = data.describe().loc[['mean', 'std'], 'edep'].values
    label = 'Layer {}'.format(layer)
#            + '\n' + r'$\mu = {:.3f}$ MeV'.format(mean)
#            + '\n' + r'$\sigma = {:.3f}$ MeV'.format(std)
    ax.text(*bpos, label, transform=ax.transAxes, fontsize=26, verticalalignment='top', bbox=bprops)
fig.text(0.5, 0.02, 'Energy Deposit (MeV)', fontsize=26, ha='center')
fig.text(0.02, 0.5, 'Normalized Counts', fontsize=26, va='center', rotation=90)
fig.savefig(os.path.join(args.outdir, 'edep_layers.png'))


if args.save:
    prof.to_csv(args.save)

