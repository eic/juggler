'''
    A utility script to help the analysis of imaging calorimeter data

    Author: Chao Peng (ANL)
    Date: 04/30/2021
'''

import os
import ROOT
import numpy as np
import pandas as pd
import matplotlib
import DDG4
from ROOT import gROOT, gInterpreter


# helper function to truncate color map (for a better view from the rainbow colormap)
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


# load root macros, input is an argument string
def load_root_macros(arg_macros):
    for path in arg_macros.split(','):
        path = path.strip()
        if os.path.exists(path):
            gROOT.Macro(path)
        else:
            print('\"{}\" does not exist, skip loading it.'.format(path))


# read mc particles from root file
def get_mcp_data(path, evnums=None, branch='mcparticles2'):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(2000*len(evnums), 6))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        # extract full mc particle data
        for part in getattr(events, branch):
            dbuf[idb] = (iev, part.psx, part.psy, part.psz, part.pdgID, part.status)
            idb += 1
    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'px', 'py', 'pz', 'pid', 'status'])


# read hits data from root file
def get_hits_data(path, evnums=None, branch='EcalBarrelClustersLayers'):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(2000*len(evnums), 7))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        for layer in getattr(events, branch):
            for k, hit in enumerate(layer.hits):
                if k < layer.nhits:
                    dbuf[idb] = (iev, layer.clusterID, layer.layerID, hit.x, hit.y, hit.z, hit.E)
                    idb += 1

    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'cluster', 'layer', 'x', 'y', 'z', 'edep'])


# read layers data from root file
def get_layers_data(path, evnums=None, branch="EcalBarrelClustersLayers"):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(2000*len(evnums), 7))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        for layer in getattr(events, branch):
            dbuf[idb] = (iev, layer.clusterID, layer.layerID, layer.position.x, layer.position.y, layer.position.z, layer.edep)
            idb += 1

    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'cluster', 'layer', 'x', 'y', 'z', 'edep'])


def compact_constants(path, names):
    if not os.path.exists(path):
        print('Cannot find compact file \"{}\".'.format(path))
        return []
    kernel = DDG4.Kernel()
    description = kernel.detectorDescription()
    kernel.loadGeometry("file:{}".format(path))
    try:
        vals = [description.constantAsDouble(n) for n in names]
    except:
        print('Fail to extract values from {}, return empty.'.format(names))
        vals = []
    kernel.terminate()
    return vals

