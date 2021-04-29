import os
import numpy as np
import pandas as pd
import ROOT
from ROOT import gROOT, gInterpreter
import argparse
import matplotlib
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


# define Ecal shapes
rmin, thickness, length = 890, 20*(10. + 1.65), 860*2+500

parser = argparse.ArgumentParser(description='sampling calorimeter cluster imaging')
parser.add_argument('file', type=str, help='path to root file')
parser.add_argument('-e', type=int, default=1, dest='iev', help='event number to check')
parser.add_argument('-o', type=str, default="./plots", dest='outdir', help='output directory')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

gROOT.Macro('rootlogon.C')
f = ROOT.TFile(args.file)
events = f.events

if args.iev >= events.GetEntries():
    print('Error: event {:d} is out of range (0 - {:d})'.format(args.iev, events.GetEntries() - 1))
    exit(-1)

# prepare data
events.GetEntry(args.iev)
data = []
for i, layer in enumerate(events.EcalBarrelClustersLayers):
    if layer.clusterID > 0:
        continue
    for k, hit in enumerate(layer.hits):
        if k < layer.nhits:
            data.append((layer.layerID, hit.x, hit.y, hit.z, hit.E))
df = pd.DataFrame(data=data, columns=['layer', 'x', 'y', 'z', 'edep'])
# print(df)

mc_p = (events.mcparticles2[2].psx, events.mcparticles2[2].psy, events.mcparticles2[2].psz)


# cluster plot
fig = plt.figure(figsize=(20, 16), dpi=160)
ax = fig.add_subplot(111, projection='3d')

# draw particle line
vec = mc_p/np.linalg.norm(mc_p)
length = rmin/np.sqrt(vec[0]**2 + vec[1]**2)
ax.plot([0., vec[2]*length], [0., vec[1]*length], '--', zs=[0., vec[0]*length], color='green')

# normalize to get colors
cmap = plt.get_cmap('rainbow')
emin = df['edep'].min()
emax = df['edep'].max()
vals = (df['edep'] - emin)/(emax - emin)
colors = cmap(vals)

# hits
ax.scatter(df['z'], df['y'], df['x'], c=colors, marker='o', s=1.0)
ax.tick_params(labelsize=24)
ax.set_xlabel('z (mm)', fontsize=26, labelpad=20)
ax.set_ylabel('y (mm)', fontsize=26, labelpad=20)
ax.set_zlabel('x (mm)', fontsize=26, labelpad=20)
cb = plt.colorbar(cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=df['edep'].min(), vmax=df['edep'].max()), cmap=cmap),
                  ax=ax, shrink=0.85)
cb.ax.tick_params(labelsize=24)
cb.ax.get_yaxis().labelpad = 15
cb.ax.set_ylabel('Energy Deposit (MeV)', rotation=90, fontsize=28)

# cylinder dimension (detection plane)
x = np.linspace(-rmin, rmin, 500)
x2 = np.linspace(-(rmin + thickness), rmin + thickness, 500)
z = np.linspace(-length, length, 500)
Xc, Zc = np.meshgrid(x, z)
Xc2, Zc2 = np.meshgrid(x2, z)
Yc = np.sqrt(rmin**2 - Xc**2)
Yc2 = np.sqrt((rmin+thickness)**2 - Xc2**2)

rstride = 20
cstride = 10
ax.plot_surface(Zc, Yc, Xc, alpha=0.1, rstride=rstride, cstride=cstride, color='royalblue')
ax.plot_surface(Zc, -Yc, Xc, alpha=0.1, rstride=rstride, cstride=cstride, color='royalblue')
ax.plot_surface(Zc2, Yc2, Xc2, alpha=0.1, rstride=rstride, cstride=cstride, color='forestgreen')
ax.plot_surface(Zc2, -Yc2, Xc2, alpha=0.1, rstride=rstride, cstride=cstride, color='forestgreen')
ax.set_zlim(-(rmin + thickness), rmin + thickness)
ax.set_xlim(-(rmin + thickness), rmin + thickness)
ax.set_ylim(-(rmin + thickness), rmin + thickness)
fig.tight_layout()
fig.savefig(os.path.join(args.outdir, 'e{}_cluster.png'.format(args.iev)))

# zoomed-in plot
fig = plt.figure(figsize=(20, 16), dpi=160)
ax = fig.add_subplot(111, projection='3d')
ax.plot([0., vec[2]*length], [0., vec[1]*length], '--', zs=[0., vec[0]*length], color='green')
ax.scatter(df['z'], df['y'], df['x'], c=colors, marker='o', s=30.0)
ax.set_xlim(vec[2]*length - 100., vec[2]*length + 100.)
ax.set_ylim(vec[1]*length - 100., vec[1]*length + 200.)
ax.set_zlim(vec[0]*length - 100., vec[0]*length + 200.)
ax.tick_params(labelsize=24)
ax.set_xlabel('z (mm)', fontsize=26, labelpad=20)
ax.set_ylabel('y (mm)', fontsize=26, labelpad=20)
ax.set_zlabel('x (mm)', fontsize=26, labelpad=20)
cb = plt.colorbar(cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=df['edep'].min(), vmax=df['edep'].max()), cmap=cmap),
                  ax=ax, shrink=0.85)
cb.ax.tick_params(labelsize=24)
cb.ax.get_yaxis().labelpad = 15
cb.ax.set_ylabel('Energy Deposit (MeV)', rotation=90, fontsize=28)

# data = df.groupby('layer')[['z', 'y','x']].mean().values
# ax.scatter(*data.T, s=50.0, color='black')
data = df[['z', 'y', 'x']].values
datamean = data.mean(axis=0)
uu, dd, vv = np.linalg.svd(data - datamean)
linepts = vv[0] * np.mgrid[-200:200:2j][:, np.newaxis]
linepts += datamean
ax.plot3D(*linepts.T, 'k:')

fig.tight_layout()
fig.savefig(os.path.join(args.outdir, 'e{}_cluster_zoom.png'.format(args.iev)))


# projection plot
df['phi'] = np.arctan2(df['y'].values, df['x'].values)
df['r'] = np.sqrt(df['x'].values**2 + df['y'].values**2 + df['z'].values**2)
df['theta'] = np.arccos(df['z'].values/df['r'].values)
vecp = (np.arccos(vec[2]), np.arctan2(vec[1], vec[0]))

tanth = np.tan(df['theta'].values)
sinph = np.sin(df['phi'].values)
cosph = np.cos(df['phi'].values)
energies = df['edep'].values

xvals = np.abs(tanth)*cosph
yvals = np.abs(tanth)*sinph

# color map
cmap = truncate_colormap(plt.get_cmap('jet'), 0.1, 0.9)

fig, axs = plt.subplots(1, 2, figsize=(17, 16), dpi=160, gridspec_kw={'wspace':0., 'width_ratios': [16, 1]})
ax = axs[0]
# cluster angle range
phi_rg = ((vecp[1] - 0.05)*1000., (vecp[1] + 0.05)*1000.)
th_rg = ((vecp[0] - 0.05)*1000., (vecp[0] + 0.05)*1000.)
# h = ax.hist2d(xvals, yvals, weights=energies, bins=(np.linspace(-0.1, 0.1, 100), np.linspace(-0.1, 0.1, 100)), cmin=df['edep'].min())
h = ax.hist2d(df['theta'].values*1000., df['phi'].values*1000., weights=energies,
              bins=(np.arange(*th_rg, step=1.), np.arange(*phi_rg, step=1.)),
              cmap=cmap, cmin=df['edep'].min())
ax.set_ylabel(r'$\phi$ (mrad)', fontsize=28)
ax.set_xlabel(r'$\theta$ (mrad)', fontsize=28)
ax.set_ylim(*phi_rg)
ax.set_xlim(*th_rg)
ax.tick_params(labelsize=24)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(10))
ax.grid(linestyle=':', which='both')
ax.set_axisbelow(True)
cb = fig.colorbar(h[3], cax=axs[1], pad=0.)
cb.ax.tick_params(labelsize=24)
cb.ax.get_yaxis().labelpad = 15
cb.ax.set_ylabel('Energy Deposit (MeV)', rotation=90, fontsize=28)
# print(df[['theta', 'phi']])
# print(xvals, yvals)
# print(df['phi'].describe())
fig.savefig(os.path.join(args.outdir, 'e{}_topo.png'.format(args.iev)))

