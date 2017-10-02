from statsmodels.sandbox.stats.multicomp import multipletests
import rpy2.robjects as robjects
import pandas as pd
from collections import Counter
import operator
import matplotlib.pyplot as plt
import numpy as np
import math
import ipywidgets as widgets
import os

def dropdown():
    wd = widgets.Dropdown(
        options={'Gene Ontology: biological process': 'go_bp.gmt', 
                 'Gene Ontology: molecular function': 'GOmf', 
                 'Gene Ontology: cellular component': 'GOcc'},
        value='go_bp.gmt',
        description='Select the database to perform the functional enrichment:',
    )
    return wd

def dropdown_gene_set(allenrichment_results):
    wd = widgets.Dropdown(
        options = {gs:gs for gs in sorted(allenrichment_results.TERM.tolist())},
        value = sorted(allenrichment_results.TERM.tolist())[0],
        description= 'Select a gene set to check the genes contained in it'
    )
    return wd

def filename_groups():
    d = {fnam:fnam for fnam in os.listdir('../data/') if os.path.isfile('../data/'+fnam) and fnam.startswith('.')==False}
    d['No_file'] = 'No_file'
    wd = widgets.Dropdown(
        options = d,
        value = 'No_file',
        description= 'Select the groups file'
    )
    return wd 

def filename_GO():
    d = {fnam:fnam for fnam in os.listdir('../data/') if os.path.isfile('../data/'+fnam) and fnam.startswith('.')==False}
    d['No_file'] = 'No_file'
    wd = widgets.Dropdown(
        options = d,
        value = 'No_file',
        description= 'Select a file with list of gene sets to enrich:'
    )
    return wd 

def pval_thres():
    fs = widgets.FloatSlider(
        value=0.05,
        min=0,
        max=0.5,
        step=0.05,
        description='Corrected P-value (Q-value) threshold:',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.2f',
    )
    return fs
        
def load_data(filename,db_name):
    #Load cluster data
    clusters = pd.read_excel('../data/'+filename)
    cluster_genes = {}
    for cluster in clusters.columns.tolist():
        cluster_genes[cluster] = clusters[cluster].dropna().tolist()
    #Load functional annotation data
    genesets = {}
    allgenes = []
    ffile = open('../data/MSigDB/'+db_name,'r')
    for line in ffile.readlines():
        f = line.strip().split('\t')
        genesets[f[0]] = f[2:]
        allgenes += f[2:]
    ffile.close()
    allgenes = set(allgenes)
    return allgenes,genesets,cluster_genes,clusters

def functional_enrichment(term_genes, N, k_l, alpha):
    functional_data = pd.DataFrame()
    k = len(k_l)
    term_data = []
    for term in term_genes.keys():
        m = len(term_genes[term])
        xl = list(set(term_genes[term]).intersection(set(k_l)))
        x = len(xl)

        if x != 0:
            xlist = []
            for i in range(x, m + 1):
                xlist.append(i)

            # calculation of the hypervalue
            dhyper = robjects.r['dhyper']
            xboh = robjects.IntVector(xlist)
            dhypervalue = dhyper(xboh, m, (N - m), k, log=False)
            # print(x,N,m,k,dhypervalue)
            term_data.append({'TERM': term, 'PVALUE': sum(dhypervalue)})

    term_data_df = pd.DataFrame(term_data)

    if len(term_data_df) > 0:
        corrected_p = multipletests(term_data_df.PVALUE.tolist(), method='fdr_bh', is_sorted=False, returnsorted=False)
        term_data_df['QVALUE'] = list(corrected_p[1])
        functional_data = functional_data.append(term_data_df)

        return functional_data[functional_data['QVALUE'] < alpha]
    else:
        return ''

def enrichment_all_groups(genesets,filt,allgenes,cluster_genes,Pvalue_thr):
    if filt != 'No_file':
        gos_sel = pd.read_excel('../data/'+filt,header=None,names=['TERM'])
        term_genes = {k:v for k,v in genesets.items() if k in gos_sel['TERM'].tolist()}
    else:
        term_genes = genesets


    allenrichment_results = pd.DataFrame()
    for cluster, cluster_genes_l in cluster_genes.items():
        print('Working on....',cluster)
        enrichment_results = functional_enrichment(term_genes=term_genes, N=len(allgenes), k_l=cluster_genes_l, alpha=Pvalue_thr)
        enrichment_results['GROUP'] = cluster
        allenrichment_results = allenrichment_results.append(enrichment_results)
    return allenrichment_results

def plot_heatmap(allenrichment_results,clusters,cluster_genes,Pvalue_thr):
    #Heatmap plotting
    terms = dict(Counter(allenrichment_results.TERM.tolist()))
    sorted_terms = [t for t,v in sorted(terms.items(), key=operator.itemgetter(1),reverse=False)]

    heatmap_l = []
    for term in sorted_terms:
        d = {}
        enrichment_term = allenrichment_results[allenrichment_results['TERM']==term]
        for cluster in cluster_genes.keys():
            if cluster in enrichment_term.GROUP.tolist():
                d[cluster] = enrichment_term[enrichment_term['GROUP']==cluster]['QVALUE'].tolist()[0]
            else:
                d[cluster] = float('nan')
        heatmap_l.append(d)
    heatmapdf = pd.DataFrame(heatmap_l)
    
    fig, ax = plt.subplots(figsize=(len(clusters.columns.tolist()), len(sorted_terms)/3))

    my_cmap_r = plt.cm.get_cmap('YlOrRd_r', 100)
    m = np.ma.masked_where(np.isnan(heatmapdf[clusters.columns.tolist()]),heatmapdf[clusters.columns.tolist()])

    heatmap = plt.pcolor(m, cmap=my_cmap_r, linewidth=1, edgecolor="white"
                         ,vmin=heatmapdf.min().min(),vmax=Pvalue_thr)
    plt.xticks([v+0.5 for v in range(0,len(clusters.columns.tolist()))],clusters.columns.tolist(),rotation=90)
    plt.yticks([v+0.5 for v in range(0,len(sorted_terms))],sorted_terms)
    ax.xaxis.set_tick_params(labeltop='on')
    plt.colorbar(label='Q-value scale')  
    plt.show()
    
    return