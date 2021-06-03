'''
    A script to visualize the cluster
    It reads the output from the Juggler component ImagingClusterReco, which is supposed to be clusters of hits after
    digitization, reconstruction, and clustering

    Author: Chao Peng (ANL)
    Date: 04/30/2021
'''

import os
import numpy as np
import pandas as pd
import argparse
import matplotlib
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from utils import *
import sys


# draw cluster in a 3d axis, expect a numpy array of (nhits, 4) shape with each row contains (x, y, z, E)
# note z and x axes are switched
def draw_hits3d(axis, data, cmap, units=('mm', 'mm', 'mm', 'MeV'), fontsize=24, **kwargs):
    # normalize energy to get colors
    x, y, z, edep = np.transpose(data)
    cvals = edep - min(edep) / (max(edep) - min(edep))
    cvals[np.isnan(cvals)] = 1.0
    colors = cmap(cvals)

    # hits
    axis.scatter(z, y, x, c=colors, marker='o', **kwargs)
    axis.tick_params(labelsize=fontsize)
    axis.set_zlabel('x ({})'.format(units[2]), fontsize=fontsize + 2, labelpad=fontsize)
    axis.set_ylabel('y ({})'.format(units[1]), fontsize=fontsize + 2, labelpad=fontsize)
    axis.set_xlabel('z ({})'.format(units[0]), fontsize=fontsize + 2, labelpad=fontsize)
    cb = plt.colorbar(cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=min(edep), vmax=max(edep)), cmap=cmap),
                      ax=axis, shrink=0.85)
    cb.ax.tick_params(labelsize=fontsize)
    cb.ax.get_yaxis().labelpad = fontsize
    cb.ax.set_ylabel('Energy Deposit ({})'.format(units[3]), rotation=90, fontsize=fontsize + 4)
    return axis


# draw a cylinder in 3d axes
# note z and x axes are switched
def draw_cylinder3d(axis, r, z, order=['x', 'y', 'z'], rsteps=500, zsteps=500, **kwargs):
    x = np.linspace(-r, r, rsteps)
    z = np.linspace(-z, z, zsteps)
    Xc, Zc = np.meshgrid(x, z)
    Yc = np.sqrt(r**2 - Xc**2)

    axis.plot_surface(Zc, Yc, Xc, alpha=0.1, **kwargs)
    axis.plot_surface(Zc, -Yc, Xc, alpha=0.1, **kwargs)
    return axis


# fit the track of cluster and draw the fit
def draw_track_fit(axis, dfh, length=200, stop_layer=8, scat_kw=dict(), line_kw=dict()):
    dfh = dfh[dfh['layer'] <= stop_layer]
    data = dfh.groupby('layer')[['z', 'y','x']].mean().values
    # data = dfh[['z', 'y', 'x']].values
    # ax.scatter(*data.T, **scat_kw)
    datamean = data.mean(axis=0)
    uu, dd, vv = np.linalg.svd(data - datamean)
    linepts = vv[0] * np.mgrid[-length:length:2j][:, np.newaxis]
    linepts += datamean
    axis.plot3D(*linepts.T, 'k:')
    return axis


# color map
def draw_heatmap(axis, x, y, weights, bins=1000, cmin=0., cmap=plt.get_cmap('rainbow'), pc_kw=dict()):
    w, xedg, yedg = np.histogram2d(x, y, weights=weights, bins=bins)
    xsz = np.mean(np.diff(xedg))
    ysz = np.mean(np.diff(yedg))
    wmin, wmax = w.min(), w.max()
    recs, clrs = [], []
    for i in np.arange(len(xedg) - 1):
        for j in np.arange(len(yedg) - 1):
            if w[i][j] > cmin:
                recs.append(Rectangle((xedg[i], yedg[j]), xsz, ysz))
                clrs.append(cmap((w[i][j] - wmin) / (wmax - wmin)))
    axis.add_collection(PatchCollection(recs, facecolor=clrs, **pc_kw))
    axis.set_xlim(xedg[0], xedg[-1])
    axis.set_ylim(yedg[0], yedg[-1])
    return axis, cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=wmin, vmax=wmax), cmap=cmap)


# execute this script
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize the cluster from analysis')
    parser.add_argument('file', type=str, help='path to root file')
    parser.add_argument('-e', type=int, default=0, dest='iev', help='event number to plot')
    parser.add_argument('-c', type=int, default=0, dest='icl', help='cluster number to plot (0: all, -1: no cluster)')
    parser.add_argument('-s', type=int, default=8, dest='stop', help='stop layer for track fit')
    parser.add_argument('-o', type=str, default='./plots', dest='outdir', help='output directory')
    parser.add_argument('--compact', type=str, default='', dest='compact', help='compact file')
    parser.add_argument('-m', '--macros', type=str, default='rootlogon.C', dest='macros',
                        help='root macros to load (accept multiple paths separated by \",\")')
    parser.add_argument('-b', '--branch-name', type=str, default='RecoEcalBarrelHits', dest='branch',
                        help='branch name in the root file (outputLayerCollection from ImagingClusterReco)')
    parser.add_argument('--topo-size', type=float, default=2.0, dest='topo_size',
                        help='bin size for projection plot (mrad)')
    parser.add_argument('--topo-range', type=float, default=50.0, dest='topo_range',
                        help='half range for projection plot (mrad)')
    parser.add_argument('--zoom-factor', type=float, default=1.0, dest='zoom_factor',
                        help='factor for zoom-in')
    args = parser.parse_args()


    # we can read these values from xml file
    desc = compact_constants(args.compact, [
        'cb_ECal_RMin',
        'cb_ECal_ReadoutLayerThickness',
        'cb_ECal_ReadoutLayerNumber',
        'cb_ECal_Length'
    ])
    if not len(desc):
        # or define Ecal shapes
        rmin, thickness, length = 890, 20*(10. + 1.65), 860*2+500
    else:
        # convert cm to mm
        rmin = desc[0]*10.
        thickness = desc[1]*desc[2]*10.
        length = desc[3]*10.


    # read data
    load_root_macros(args.macros)
    df = get_hits_data(args.file, args.iev, branch=args.branch)
    if args.icl != 0:
        df = df[df['cluster'] == args.icl]
    if not len(df):
        print("Error: do not find any hits for cluster {:d} in event {:d}".format(args.icl, args.iev))
        exit(-1)
    # convert to polar coordinates (mrad), and stack all r values
    df['r'] = np.sqrt(df['x'].values**2 + df['y'].values**2 + df['z'].values**2)
    df['phi'] = np.arctan2(df['y'].values, df['x'].values)*1000.
    df['theta'] = np.arccos(df['z'].values/df['r'].values)*1000.
    df['eta'] = -np.log(np.tan(df['theta'].values/1000./2.))

    # truth
    dfmcp = get_mcp_simple(args.file, args.iev, 'mcparticles2').iloc[0]
    pdgbase = ROOT.TDatabasePDG()
    inpart = pdgbase.GetParticle(int(dfmcp['pid']))
    print("Incoming particle = {}, pdgcode = {}, charge = {}, mass = {}"\
          .format(inpart.GetName(), inpart.PdgCode(), inpart.Charge(), inpart.Mass()))
    # neutral particle, no need to consider magnetic field
    if np.isclose(inpart.Charge(), 0., rtol=1e-5):
        vec = dfmcp[['px', 'py', 'pz']].values
    # charge particle, use the cluster center
    else:
        flayer = df[df['layer'] == df['layer'].min()]
        vec = flayer[['x', 'y', 'z']].mean().values
    vec = vec/np.linalg.norm(vec)

    # particle line from (0, 0, 0) to the inner Ecal surface
    length = rmin/np.sqrt(vec[0]**2 + vec[1]**2)
    pline = np.transpose(vec*np.mgrid[0:length:2j][:, np.newaxis])
    cmap = truncate_colormap(plt.get_cmap('jet'), 0.1, 0.9)

    os.makedirs(args.outdir, exist_ok=True)
    # cluster plot
    fig = plt.figure(figsize=(15, 12), dpi=160)
    ax = fig.add_subplot(111, projection='3d')
    # draw particle line
    ax.plot(*pline[[2, 1]], '--', zs=pline[0], color='green')
    # draw hits
    draw_hits3d(ax, df[['x', 'y', 'z', 'edep']].values, cmap, s=5.0)
    draw_cylinder3d(ax, rmin, length, rstride=10, cstride=10, color='royalblue')
    draw_cylinder3d(ax, rmin + thickness, length, rstride=10, cstride=10, color='forestgreen')
    ax.set_zlim(-(rmin + thickness), rmin + thickness)
    ax.set_ylim(-(rmin + thickness), rmin + thickness)
    ax.set_xlim(-length, length)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, 'e{}_cluster.png'.format(args.iev)))


    # zoomed-in cluster plot
    fig = plt.figure(figsize=(15, 12), dpi=160)
    ax = fig.add_subplot(111, projection='3d')
    # draw particle line
    ax.plot(*pline[[2, 1]], '--', zs=pline[0], color='green')
    # draw hits
    draw_hits3d(ax, df[['x', 'y', 'z', 'edep']].values, cmap, s=20.0)
    # draw_track_fit(ax, df, stop_layer=args.stop,
    #                scat_kw=dict(color='k', s=50.0), line_kw=dict(linestyle=':', color='k', lw=3))
    # view range
    center = (length + thickness/2.)*vec
    ranges = np.vstack([center - thickness/args.zoom_factor, center + thickness/args.zoom_factor]).T
    ax.set_zlim(*ranges[0])
    ax.set_ylim(*ranges[1])
    ax.set_xlim(*ranges[2])

    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, 'e{}_cluster_zoom.png'.format(args.iev)))


    # projection plot
    # convert to mrad
    vecp = np.asarray([np.arccos(vec[2]), np.arctan2(vec[1], vec[0])])*1000.
    phi_rg = np.asarray([vecp[1] - args.topo_range, vecp[1] + args.topo_range])
    th_rg = np.asarray([vecp[0] - args.topo_range, vecp[0] + args.topo_range])
    eta_rg = np.resize(-np.log(np.tan(vecp[0]/1000./2.)), 2) + np.asarray([-args.topo_range, args.topo_range])/1000.

    fig, axs = plt.subplots(1, 2, figsize=(13, 12), dpi=160, gridspec_kw={'wspace':0., 'width_ratios': [12, 1]})
    ax, sm = draw_heatmap(axs[0], df['eta'].values, df['phi'].values, weights=df['edep'].values,
                          bins=(np.arange(*eta_rg, step=args.topo_size/1000.), np.arange(*phi_rg, step=args.topo_size)),
                          cmap=cmap, cmin=0., pc_kw=dict(alpha=0.8, edgecolor='k'))

    ax.set_ylabel(r'$\phi$ (mrad)', fontsize=32)
    ax.set_xlabel(r'$\eta$', fontsize=32)
    ax.tick_params(labelsize=28)
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.grid(linestyle=':', which='both')
    ax.set_axisbelow(True)
    cb = plt.colorbar(sm, cax=axs[1], shrink=0.85, aspect=1.2*20)
    cb.ax.tick_params(labelsize=28)
    cb.ax.get_yaxis().labelpad = 10
    cb.ax.set_ylabel('Energy Deposit (MeV)', rotation=90, fontsize=32)
    fig.savefig(os.path.join(args.outdir, 'e{}_topo.png'.format(args.iev)))

