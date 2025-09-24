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

def get_gene_data(data, gene_id): #Dato l'AnnData e l'id del gene, ritorna quante volte compare quel gene PER ogni spot
    array_data = data.X.toarray().T
    gene_data  = array_data[gene_id]
    return gene_data

def select_top_spatially_variable_genes(data, n_top, min_gene_expression=300, n_top_genes=3000): #La funzione mi ritorna gli indici dei geni più variabili
    # n_top = quanti geni vuoi tenere alla fine.

    # min_gene_expression = soglia minima di espressione totale (somma su tutte le celle/spot) per non scartare un gene.

    # n_top_genes = quanti geni far considerare all’algoritmo HVG di Scanpy (Seurat v3) nella fase di ranking.


    tot_gene_expression = data.X.toarray().sum(axis=0) #Calcolo per ogni gene quante volte compare in totale (non per spot, ma in totale)
    
    # sc.pp.normalize_total(data, target_sum=1e4)
    # sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=n_top_genes) #Fa l'HVG e lo mette in Data. la lista è dal migliore (più variabile) al peggiore

    var_rank = data.var["highly_variable_rank"].to_numpy()
    var_rank = np.nan_to_num(var_rank, nan=np.nanmax(var_rank)+1) 
    max_rank = np.nanmax(var_rank) #Peggior rank possibile

    low_expression_indexes = np.where(tot_gene_expression < min_gene_expression)[0] #Prendo gli indici dei geni che hanno un'espressione totale minore di 300
    var_rank[low_expression_indexes] = max_rank #A questi geni gli do il peggior rank possibile, così non li considero

    smallest = np.argsort(var_rank)[:n_top] #Ordina gli indici dei geni per rank crescente (quindi prima i più variabili).
    smallest = smallest[(var_rank[smallest] != max_rank)]

    return smallest #array NumPy di indici di colonna (riferiti a data.var/data.X) dei geni selezionati come “top variabili”