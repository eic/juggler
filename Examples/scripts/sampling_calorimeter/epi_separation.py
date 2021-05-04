'''
    A script to use the energy profile for e-pi separation
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


def prepare_data(path, **kwargs):
    tmp = get_layers_data(path, **kwargs)
    tmp.loc[:, 'total_edep'] = tmp.groupby(['event', 'cluster'])['edep'].transform('sum')
    tmp.loc[:, 'efrac'] = tmp['edep']/tmp['total_edep']
    # tmp = tmp[tmp['cluster'] == 0]
    return tmp


def calc_chi2(grp, pr, lmin=5, lmax=12):
    grp2 = grp[(grp['layer'] >= lmin) & (grp['layer'] <= lmax)]
    mean, std = pr.loc[grp2['layer'], ['mean', 'std']].values.T*pr['energy'].mean()
    edeps = grp2['edep'].values
    return np.sqrt(np.sum((edeps - mean)**2/std**2)/float(len(edeps)))


# execute this script
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='epi_separation')
    parser.add_argument('efile', type=str, help='path to root file (electrons)')
    parser.add_argument('pifile', type=str, help='path to root file (pions)')
    parser.add_argument('--prof', type=str, default='profile.csv', help='path to electron profile')
    parser.add_argument('--plot-dir', type=str, default='./plots', dest='outdir', help='output directory')
    parser.add_argument('-b', '--branch-name', type=str, default='EcalBarrelClustersLayers', dest='branch',
                        help='branch name in the root file (outputLayerCollection from ImagingClusterReco)')
    parser.add_argument('-m', '--macros', type=str, default='rootlogon.C', dest='macros',
                        help='root macros to load (accept multiple paths separated by \",\")')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)


    load_root_macros(args.macros)
    # prepare data
    dfe = prepare_data(args.efile, branch=args.branch)
    dfpi = prepare_data(args.pifile, branch=args.branch)

    colors = ['royalblue', 'indianred', 'limegreen']
    prof = pd.read_csv(args.prof).set_index('layer', drop=True)

    # profile comparison
    fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
    for title, color in zip(sorted(prof['type'].unique()), colors):
        pr = prof[prof['type'] == title]
        ax.plot(pr.index, pr['mean'], color=color, lw=2)
        ax.fill_between(pr.index, (pr['mean'] - pr['std']), (pr['mean'] + pr['std']),
                        color=color, alpha=0.3, label=title)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.grid(linestyle=':', which='both')
    ax.tick_params(labelsize=24)
    ax.set_xlabel('Layer', fontsize=26)
    ax.set_ylabel('Energy Deposit Percentage', fontsize=26)
    ax.legend(fontsize=26, loc='upper left')
    fig.savefig(os.path.join(args.outdir, 'compare_prof.png'))


    # check profile
    layer_range = (4, 12)
    pre = prof[prof['type'].str.lower() == 'em']
    fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
    chi2 = dfe.groupby(['event', 'cluster']).apply(lambda x: calc_chi2(x, pre, *layer_range))
    ax.hist(chi2, bins=np.linspace(0, 5, 200), weights=[1/float(len(chi2))]*chi2,
            ec='royalblue', color='royalblue', alpha=0.5, label='electrons')
    chi2 = dfpi.groupby(['event', 'cluster']).apply(lambda x: calc_chi2(x, pre, *layer_range))
    # print(chi2[chi2 < 0.7])
    ax.hist(chi2, bins=np.linspace(0, 5, 200), weights=[1/float(len(chi2))]*chi2,
            ec='indianred', color='indianred', alpha=0.5, label='pions')
    ax.grid(linestyle=':', which='major')
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=24)
    ax.set_xlabel(r'$\chi^2$', fontsize=26)
    ax.set_ylabel('Normalized Counts', fontsize=26)
    # ax.set_yscale('log')
    # ax.set_ylim(1e-3, 1.)
    ax.legend(title=r'$\chi^2$ for $E_{{dep}}$ in layer {:d} - {:d}'.format(*layer_range), title_fontsize=28, fontsize=26)
    fig.savefig(os.path.join(args.outdir, 'efrac_chi2.png'))

