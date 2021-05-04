'''
    A script to generate the energy profile (layer-wise)
    It reads the output from the Juggler component ImagingClusterReco, which is supposed to be clusters of hits after
    digitization, reconstruction, and clustering

    Author: Chao Peng (ANL)
    Date: 04/30/2021
'''

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
from utils import *


def find_start_layer(grp, min_edep=0.5):
    ids, edeps = grp.sort_values('layer')[['layer', 'edep']].values.T
    edeps = np.cumsum(edeps)
    return min(ids[edeps > min_edep]) if sum(edeps > min_edep) > 0 else 20


# execute this script
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate energy profiles')
    parser.add_argument('file', type=str, help='path to root file')
    parser.add_argument('--plot-dir', type=str, default='./plots', dest='outdir', help='output directory')
    parser.add_argument('--type', type=str, default='unknown', dest='type', help='profile type (used in save)')
    parser.add_argument('--energy', type=float, default=5000., dest='energy', help='incident particle energy (MeV)')
    parser.add_argument('--save', type=str, default='', dest='save', help='path to save profile')
    parser.add_argument('--color', type=str, default='royalblue', dest='color', help='colors for bar plots')
    parser.add_argument('-b', '--branch-name', type=str, default='EcalBarrelClustersLayers', dest='branch',
                        help='branch name in the root file (outputLayerCollection from ImagingClusterReco)')
    parser.add_argument('-m', '--macros', type=str, default='rootlogon.C', dest='macros',
                        help='root macros to load (accept multiple paths separated by \",\")')
    args = parser.parse_args()

    load_root_macros(args.macros)
    # prepare data
    dfe = get_layers_data(args.file, branch=args.branch)
    dfe.loc[:, 'total_edep'] = dfe.groupby(['event', 'cluster'])['edep'].transform('sum')
    # dfe.loc[:, 'efrac'] = dfe['edep']/dfe['total_edep']
    dfe.loc[:, 'efrac'] = dfe['edep']/args.energy
    dfe = dfe[dfe['cluster'] == 0]

    os.makedirs(args.outdir, exist_ok=True)

    slayer = dfe.groupby('event').apply(lambda x: find_start_layer(x, 1.0)).astype(int)
    dfe = dfe.merge(slayer.to_frame(name='slayer'), on='event')
    # dfe.loc[:, 'eff_layer'] = dfe['layer'] - dfe['slayer']
    # prof = dfe[dfe['eff_layer'] > 0].groupby('eff_layer')['edep'].describe()
    prof = dfe.groupby('layer')['efrac'].describe()

    # print(prof['mean'])
    # plot profiles
    bpos, bprops = (0.5, 0.95), dict(boxstyle='round', facecolor='wheat', alpha=0.3)
    fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
    nev = len(dfe['event'].unique())
    ax.hist(dfe.groupby('event')['slayer'].min(), weights=[1/float(nev)]*nev,
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
    ax.plot(prof.index, prof['mean'].values*100., color=args.color, lw=2)
    # ax.fill_between(prof.index, prof['25%'], prof['75%'], color=args.color, alpha=0.3)
    ax.fill_between(prof.index, (prof['mean'] - prof['std'])*100., (prof['mean'] + prof['std'])*100.,
                    color=args.color, alpha=0.3)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.grid(linestyle=':', which='both')
    ax.tick_params(labelsize=24)
    ax.set_xlabel('Layer', fontsize=26)
    ax.set_ylabel('Energy Deposit Percentage', fontsize=26)
    fig.savefig(os.path.join(args.outdir, 'efrac.png'))

    layers = np.asarray([
        [1, 5, 8,],
        [10, 15, 20],
    ])

    fig, ax = plt.subplots(*layers.shape, figsize=(16, 9), dpi=160, sharex='col', sharey='all',
                           gridspec_kw=dict(hspace=0.05, wspace=0.05))

    for ax, layer in zip(ax.flat, layers.flatten()):
        data = dfe[dfe['layer'] == layer]
        ax.hist(data['efrac'].values*100., weights=[1/float(len(data))]*len(data), bins=np.linspace(0, 30, 60),
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
        ax.set_axisbelow(True)
        ax.set_yscale('log')
    fig.text(0.5, 0.02, 'Energy Deposit Percentage', fontsize=26, ha='center')
    fig.text(0.02, 0.5, 'Normalized Counts', fontsize=26, va='center', rotation=90)
    fig.savefig(os.path.join(args.outdir, 'efrac_layers.png'))


    if args.save:
        prof.loc[:, 'energy'] = args.energy
        prof.loc[:, 'type'] = args.type
        if os.path.exists(args.save):
            prev = pd.read_csv(args.save).set_index('layer', drop=True)
            prof = pd.concat([prof, prev])
        prof.to_csv(args.save)

