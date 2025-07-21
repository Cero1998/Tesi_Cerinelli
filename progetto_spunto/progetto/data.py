import scanpy as sc
import numpy  as np
from   filt import *

def load_data():
    h5ad_file       = '1.DLPFC/151673/scRNA.h5ad'
    expression_file = '1.DLPFC/151673/filtered_feature_bc_matrix.h5'
    oprefix         = 'DLPFC-151673'
    
    adata_gene_expression = sc.read_visium('1.DLPFC/151673', count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata_gene_expression.var_names_make_unique()

    init_filter_support_data(adata_gene_expression)
    
    return adata_gene_expression

def get_gene_data(data, gene_id):
    array_data = data.X.toarray().T
    gene_data  = array_data[gene_id]
    return gene_data

def select_top_spatially_variable_genes(data, n_top, min_gene_expression=300, n_top_genes=3000):

    tot_gene_expression = data.X.toarray().sum(axis=0)
    
    # sc.pp.normalize_total(data, target_sum=1e4)
    # sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=n_top_genes)

    var_rank = data.var["highly_variable_rank"].to_numpy()
    var_rank = np.nan_to_num(var_rank, nan=np.nanmax(var_rank)+1)
    max_rank = np.nanmax(var_rank)

    low_expression_indexes = np.where(tot_gene_expression < min_gene_expression)[0]
    var_rank[low_expression_indexes] = max_rank

    smallest = np.argsort(var_rank)[:n_top]
    smallest = smallest[(var_rank[smallest] != max_rank)]

    return smallest