#!/usr/bin/env python

import sys

import scanpy as sc
import numpy
import matplotlib.pyplot as plt


def main():
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/',
                            var_names='gene_symbols', cache=True)

    adata.var_names_make_unique()

    sc.tl.pca(adata, svd_solver='arpack')
    raw = adata.copy()
    #sc.pl.pca(adata, title='Unfiltered', save="_unfiltered.pdf")

    print("# cells, # genes before filtering:", adata.shape)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print("# cells, # genes after filtering:", adata.shape)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                               log1p=False, inplace=True)
    #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    #         jitter=0.4, multi_panel=True)
    
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    print("# cells, # genes after MT filtering:", adata.shape)
    #sc.pl.highest_expr_genes(adata, n_top=20)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3,
                                min_disp=0.5)
    #sc.pl.highly_variable_genes(adata)
    adata.write("filtered_data.h5")

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    print("# cells, # genes after variability filtering:", adata.shape)

    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, svd_solver='arpack')
    #sc.pl.pca_variance_ratio(adata, log=True)
    #sc.pl.pca(adata, color='CST3')

    #fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    #sc.pl.pca(raw, ax=ax[0], title="Uniltered", show=False)
    #sc.pl.pca(adata, ax=ax[1], title="Filtered", show=False)
    #plt.tight_layout()
    #plt.savefig("pca.pdf")
    #plt.close()

    adata.write('variable_data.h5')

main()
