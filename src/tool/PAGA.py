"paga_script.py"

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy as sc
data_path = "./data/"

def transpose_start_terminal_gene_exp(gene_exp_file, transpose_target_file):
    gene_exp = pd.read_table(data_path + gene_exp_file)
    gene_exp = gene_exp.transpose()
    gene_exp.to_csv(transpose_target_file, header=False)
    

def get_cell_types(annotation_file):
    annotations = pd.read_table(data_path + annotation_file)
    annotations.transpose()
    return annotations.iloc[0,:]

def generate_anndata(annotation_file, transpose_file):
    adata = sc.read_csv(transpose_file)
    adata.obs['CellType'] = get_cell_types(annotation_file)
    return adata


def paga(*args):

    start = args[0]
    data = args[2]
    annotation_file =  args[3]
    outdir = args[4]

    ## sc parameters config
    
    sc.set_figure_params(scanpy=True, frameon=False, figsize =(6,6),vector_friendly=True, fontsize=18, color_map=None, format='png', facecolor=None, transparent=False, ipython_format='png2x')
    sc.settings.figdir = outdir + "/PAGA/"

    TRANSPOSED_FILE_PATH = outdir + "/transposed_data.csv"

    transpose_start_terminal_gene_exp(data, TRANSPOSED_FILE_PATH)
    adata = generate_anndata(annotation_file, TRANSPOSED_FILE_PATH)

    sc.pp.normalize_total(adata)
    sc.tl.pca(adata, svd_solver='arpack')#,n_comps = 50)
    #sc.pl.pca_variance_ratio(adata, show = False, save = "PCA_var.png")

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=11)
    sc.tl.draw_graph(adata)
    sc.pl.draw_graph(adata,  color='CellType',legend_loc='on data',show = False, save = 'PCA.png')


    sc.tl.diffmap(adata, n_comps= 15)
    sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_diffmap')
    sc.tl.draw_graph(adata)
    sc.pl.draw_graph(adata, color='CellType', legend_loc='on data',show = False, save = 'Diffusion_map.png')

    sc.tl.louvain(adata, resolution=1.0)
    adata.obs['Louvain']= adata.obs['louvain']

    sc.tl.paga(adata, groups='Louvain')
    sc.pl.paga(adata, color=['CellType','Louvain'], save = "Louvain_clu.png", show = False)

    adata.obs['louvain'].cat.categories
    adata.obs['louvain_anno'] = adata.obs['louvain']
    sc.tl.paga(adata, groups='louvain_anno')

    sc.tl.draw_graph(adata, init_pos='paga')
    sc.pl.draw_graph(adata, color=['CellType'], legend_loc='on data', save = "Clusters.png",show = False)
    
    print("The starting cells are: ", start)
    adata.uns['iroot'] = np.flatnonzero(adata.obs['CellType']  == start)[0]


    sc.tl.dpt(adata)

    adata.obs['DPT_Pseudotime'] = adata.obs['dpt_pseudotime']

    adata.obs['b)'] = adata.obs['dpt_pseudotime']
    adata.obs['a)'] = adata.obs['CellType']

    sc.pl.paga(adata, color=['CellType', 'DPT_Pseudotime'], show = False, save = 'Pseudotime.png')

    sc.pp.log1p(adata)
    sc.pp.scale(adata)

   
    
    core = pd.read_table(outdir+"/TransSynW/cores.tsv")
    core['Gene'] = core['Gene'].str.upper()
    sg = core[core['core']=='specific']['Gene'].to_list()
    pg = core[core['core']=='pioneer']['Gene'].to_list()

    sc.pl.draw_graph(adata, color=['CellType']+sg, legend_loc='on data', save = 'Specific_genes.png', show = False)
    sc.pl.draw_graph(adata, color=['CellType']+pg, legend_loc='on data', save = "Pioneer_genes.png", show = False)

    marker = pd.read_table(outdir+"/TransSynW/markers.tsv")
    mg = marker['Gene'].to_list()

    sc.pl.draw_graph(adata, color=['CellType']+ mg, legend_loc='on data', save ="Marker_genes.png",show = False)

    sc.pl.paga_compare(
        adata, threshold=0.03, title='', right_margin=0.2, size=50, edge_width_scale=0.5,
        legend_fontsize=12, fontsize=12, frameon=False, edges=True, show = False, save = "Comparision_1.png")

    sc.pl.paga(adata, color=['CellType','louvain','dpt_pseudotime'], show = False, save = "Comparision_2.png")


    print("PAGA analysis done!!!")
 
