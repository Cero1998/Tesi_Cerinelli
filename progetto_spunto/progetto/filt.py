import numpy as np
from   sklearn.mixture import GaussianMixture

ITNI  = None
VNI   = None
NN    = None
NM    = None
ITSNI = None
HNN   = None

MEAN_FILTER_IT = 6

def init_filter_support_data(data):
    global ITNI, VNI, NN, NM, ITSNI, HNN
    neighbours = np.array([(0, 0), (0, -2), (0, +2), (-1, -1), (-1, +1), (+1, -1), (+1, +1)], dtype=int)
    
    max_col = data.obs['array_col'].max()
    max_row = data.obs['array_row'].max()
    
    coord_to_idx = np.zeros((max_row+2, max_col+2), dtype=int)-1
    
    for i, ind in enumerate(data.obs.index):
        coord_to_idx[data.obs['array_row'][ind]-1, data.obs['array_col'][ind]-1] = i

    nspot, = data.obs.index.shape
    ITNI   = np.zeros((nspot, 7), dtype=int)
    
    for i, ind in enumerate(data.obs.index):
        col = data.obs['array_col'][ind]-1
        row = data.obs['array_row'][ind]-1

        coord = np.array([row, col], dtype=int)
        
        neigh_coord = coord + neighbours
        
        ITNI[i] = coord_to_idx[neigh_coord[:,0], neigh_coord[:,1]]
    
    VNI   = (ITNI!=-1).astype(int)
    NN    = VNI.sum(axis=1)
    NM    = (VNI.T / NN.T).T
    ITSNI = ITNI[:,1:]
    HNN   = NN//2

## MAX FILTER

def max_filter(gene_data):
    global ITNI, VNI
    return np.max(gene_data[ITNI]*VNI, axis=1)

## MEAN FILTER

def super_opt_mean_filter(gene_data):
    global ITNI, NM
    return (gene_data[ITNI] * NM).sum(axis=1)

def opt_mean_filter_iterated(gene_data, it):
    for _ in range(it):
       gene_data = super_opt_mean_filter(gene_data)
    return gene_data

## EXPANSION FILTER

def expansion_filter(gene_data):
    global ITSNI, HNN
    active_neigh = gene_data[ITSNI].sum(axis=1)
    return (active_neigh >= 3).astype(int)*gene_data + (active_neigh > 3).astype(int) * (1-gene_data)

def expansion_filter_iterated(gene_data, it):
    for _ in range(it):
       gene_data = expansion_filter(gene_data)
    return gene_data

## MIXTURE of GAUSSIAN

def transform_MoG(gene_data, mean_filter_it=MEAN_FILTER_IT, exp_filter_it=3):
    
    gene_data = opt_mean_filter_iterated(gene_data, mean_filter_it)

    gm        = GaussianMixture(n_components=2, random_state=0).fit(gene_data.reshape(-1, 1))
    gene_data = gm.predict(gene_data.reshape(-1, 1))
    m1, m2    = gm.means_
    if m2 < m1:
        gene_data = 1 - gene_data

    gene_data = expansion_filter_iterated(gene_data, exp_filter_it)

    return gene_data

def apply_MoG_to_gene_idxs(data, gene_idxs, verbose=False):

    array_data = data.X.toarray().T
    
    ROWS, = gene_idxs.shape
    COLS, = array_data[0].shape
    
    transf_data = np.zeros((ROWS, COLS), dtype=int)
    
    for idx, gene_id in enumerate(gene_idxs):
        if verbose and idx % 100 == 0:
            print(f"{idx} / {ROWS}")
        gene_data        = array_data[gene_id]
        transf_data[idx] = transform_MoG(gene_data)

    return transf_data