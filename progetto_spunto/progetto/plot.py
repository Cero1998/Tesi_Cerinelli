import copy
import random
from random import choice

from   scipy.stats import norm
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from   matplotlib.colors import LinearSegmentedColormap
import numpy as np
import scipy
import scipy.sparse
import sklearn
from   sklearn.decomposition import PCA
import statistics
from   sklearn.mixture import GaussianMixture
from sklearn import metrics
from sklearn.cluster import DBSCAN

import scanpy as sc
import squidpy as sq

from filt import *
from filt import MEAN_FILTER_IT

W = '#f8f5ff'
B = '#1e1e1e'
mpl.rcParams['text.color']       = W
mpl.rcParams['axes.labelcolor']  = W
mpl.rcParams['axes.edgecolor']   = W
mpl.rcParams['axes.facecolor']   = B
mpl.rcParams['figure.facecolor'] = B
mpl.rcParams['xtick.color']      = W
mpl.rcParams['ytick.color']      = W
mpl.rcParams['font.family']      = "monospace"

colors     = [(0.12, 0, 0), (1, 1, 1)]
cmap_name  = "black_white"
BWCMAP     = LinearSegmentedColormap.from_list(cmap_name, colors)

SHAPE     = "hex"
CMAP      = "bone"
FIGSIZE   = (5, 5)
IMG_ALPHA = 0
COLORBAR  = False

def plot_spatial_scatter(
        data, 
        gene_data, 
        ax           = None, 
        cmap         = CMAP, 
        figsize      = FIGSIZE, 
        spines_color = B,
        size         = None, 
        shape        = SHAPE):
    if not ax: 
        ax = plt.gca()
    
    data.obs["show_gene"] = gene_data
    sq.pl.spatial_scatter(
            data, 
            shape     = shape,
            color     = "show_gene", 
            figsize   = figsize,
            cmap      = cmap,
            colorbar  = COLORBAR,
            img_alpha = IMG_ALPHA,
            ax        = ax,
            size      = size)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_title("")

    ax.spines['bottom'].set_color(spines_color)
    ax.spines['top'].set_color(spines_color)
    ax.spines['right'].set_color(spines_color)
    ax.spines['left'].set_color(spines_color)

def plot_transformed_MoG(data, gene_id, ax, cmap=CMAP, figsize=FIGSIZE): 
    array_data = data.X.toarray().T
    gene_data  = array_data[gene_id].copy()

    gene_data  = transform_MoG(gene_data)
    
    plot_spatial_scatter(data, gene_data, ax, cmap=cmap, figsize=figsize)

    return gene_data

def full_plot(data, gene_id, save=False, path=None):
    array_data = data.X.toarray().T
    gene_data  = array_data[gene_id]
    gene_name  = data.var.iloc[[gene_id]].index[0]
    title      = f"id: {gene_id}, name: {gene_name}"
    
    # fig, axs = plt.subplots(ncols=5, nrows=2, figsize=(25,10), facecolor=B, dpi=200)
    # ax11, ax12, ax13, ax14, ax15 = axs[0]
    # ax21, ax22, ax23, ax24, ax25 = axs[1]
    
    fig = plt.figure(figsize=(25,10), facecolor=B, dpi=100)
    
    gs  = fig.add_gridspec(2,5)
    
    ax11 = fig.add_subplot(gs[0, 0])
    ax12 = fig.add_subplot(gs[0, 1])
    ax13 = fig.add_subplot(gs[0, 2])
    
    ax21 = fig.add_subplot(gs[1, 0])
    ax22 = fig.add_subplot(gs[1, 1])
    ax23 = fig.add_subplot(gs[1, 2])

    aximg = fig.add_subplot(gs[:,3:])

    fig.suptitle(title, fontsize=20)

    #### ax12: SPATIAL SCATTER

    plot_spatial_scatter(data, gene_data, ax12)

    #### ax22: MEAN FILTERED SCATTER

    filtered_gene_data = opt_mean_filter_iterated(gene_data.copy(), MEAN_FILTER_IT)

    plot_spatial_scatter(data, filtered_gene_data, ax22)

    #### ax11: DISTRIBUTION
    
    gmax = int(gene_data.max())
        
    near_zero_idx      = gene_data < 0.4
    near_zero_perc     = near_zero_idx.astype(int).sum() / len(gene_data)
    not_zero_idx       = np.logical_not(near_zero_idx)
    gene_data_not_zero = gene_data[not_zero_idx]
    
    ax11.hist(gene_data_not_zero, bins="auto", density=True)
    ax11.set_title(f"near 0 perc.: {near_zero_perc*100:.2f}%")
    ax11.spines['top'].set_color(B)
    ax11.spines['right'].set_color(B)
    
    #### ax21: ORDERED CDF

    ax21.scatter(range(len(gene_data_not_zero)), np.sort(gene_data_not_zero), s=10, alpha=0.05)
    ax21.set_ylim([0, gmax+1])
    ax21.spines['top'].set_color(B)
    ax21.spines['right'].set_color(B)

    #### ax13 SMOOTH DISTRIBUTION

    filtered_gene_data = opt_mean_filter_iterated(gene_data.copy(), MEAN_FILTER_IT)

    ax13.hist(filtered_gene_data, bins="auto", density=True)
    ax13.spines['top'].set_color(B)
    ax13.spines['right'].set_color(B)
    
    gm = GaussianMixture(n_components=2, random_state=0).fit(filtered_gene_data.reshape(-1, 1))
    
    m1, m2 = gm.means_
    w1, w2 = gm.weights_
    s1, s2 = np.sqrt(gm.covariances_)
    
    if m1 > m2:
        tm, ts, tw = m1, s1, w1
        m1, s1, w1 = m2, s2, w2
        m2, s2, w2 = tm, ts, tw
    
    gmax = int(filtered_gene_data.max())

    ax13.axvline(m1+2.5*s1, color="r", alpha=0.5)
    ax13.axvline(m1-2.5*s1, color="r", alpha=0.5)
    
    minx = filtered_gene_data.min()
    maxx = filtered_gene_data.max()
    x    = np.linspace(minx-0.5, maxx + 0.5, 200)
    g1y  = norm(m1, s1).pdf(x)*w1
    g2y  = norm(m2, s2).pdf(x)*w2
    
    ax13.plot(x, g1y.reshape(-1), "r", alpha=0.5)
    ax13.plot(x, g2y.reshape(-1), "y", alpha=0.5)


    #### ax23: SMOOTH ORDERED CDF

    ax23.scatter(range(len(filtered_gene_data)), np.sort(filtered_gene_data), s=10, alpha=0.05)
    ax23.set_ylim([0, gmax+1])
    ax23.spines['top'].set_color(B)
    ax23.spines['right'].set_color(B)
    
    #### ax14, ax15, ax24, ax25

    # ax14.set_visible(False)
    # ax15.set_visible(False)
    # ax24.set_visible(False)
    
    # ax25.set_visible(False)
    plot_transformed_MoG(data, gene_id, aximg, cmap=BWCMAP, figsize=(10, 10))

    #### 
    if save and path:
        plt.savefig(path, format="png")
        plt.close()
    else:
        plt.show()