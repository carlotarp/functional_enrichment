from statsmodels.sandbox.stats.multicomp import multipletests
import pandas as pd
from collections import Counter
import operator
import matplotlib.pyplot as plt
import numpy as np
import math
import ipywidgets as widgets
from ipywidgets import Layout
import os
from scipy.stats import fisher_exact
import gseapy as gp
from scipy.stats import rankdata
from itertools import combinations

#GENERATE HTML PARAMETERS: first cell
def dropdown():
    print('Select the database to perform the functional enrichment:')
    wd = widgets.Dropdown(
        options={'Gene Ontology: biological process': 'go_bp.gmt', 
                 'Gene Ontology: molecular function': 'go_mf.gmt', 
                 'Gene Ontology: cellular component': 'go_cc.gmt',
                 'KEGG,Reactome & Biocarta': 'Curated_canonical_pathways_MSigDB.gmt',
                 'MSigDB hallmarks':'hallmarks.gmt',
                 'MSigDB oncogenic signatures':'oncogenic_signatures.gmt',
                 'MSigDB transcription factor targets':'transcription_factors.gmt'},
        value='go_bp.gmt',
        description='',
        layout=Layout(width='100%', height='100%')
    )
    return wd

def dropdown_gene_set(allenrichment_results):
    print('Select a gene set to check the genes contained in it')
    wd = widgets.Dropdown(
        options = {gs:gs for gs in sorted(allenrichment_results.TERM.tolist())},
        value = sorted(allenrichment_results.TERM.tolist())[0],
        description= '',
        layout=Layout(width='100%', height='100%')
    )
    return wd

def filename_groups():
    d = {fnam:fnam for fnam in os.listdir('../data/') if os.path.isfile('../data/'+fnam) and fnam.startswith('.')==False}
    d['No_file'] = 'No_file'
    print('Select the groups file')
    wd = widgets.Dropdown(
        options = d,
        value = 'No_file',
        description= '',
        layout=Layout(width='100%', height='100%')
    )
    return wd 

def filename_GO():
    d = {fnam:fnam for fnam in os.listdir('../data/') if os.path.isfile('../data/'+fnam) and fnam.startswith('.')==False}
    d['No_file'] = 'No_file'
    print('OPTIONAL - Select a file with list of gene sets to enrich (filter):')
    wd = widgets.Dropdown(
        options = d,
        value = 'No_file',
        description= '',
        layout=Layout(width='100%', height='100%')
    )
    return wd 

def min_gene_intersected():
    print('Minimum of genes intersecting gene set')
    fs = widgets.FloatSlider(
        value=0,
        min=0,
        max=100,
        step=5,
        description='',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.0f',
        layout=Layout(width='50%', height='100%'),
        position='right'
    )
    return fs

def pval_thres():
    print('Corrected P-value (Q-value) threshold')
    pv = widgets.FloatText(
        value=0.05,
        description='',
        disabled=False
    )
    return pv
################################################################################

####################### Functional enrichment functions ################################
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

def functional_enrichment(term_genes, N, k_l, alpha, xmin):
    k = len(k_l)
    term_data = []
    for term in term_genes.keys():
        m = len(term_genes[term])
        xl = list(set(term_genes[term]).intersection(set(k_l)))
        x = len(xl)
        if x >= xmin and x !=0:
            interm = [x,m-x]
            outterm = [k-x,N-m]
            odds,pvalf = fisher_exact([interm,outterm])
            term_data.append({'TERM': term, 'PVALUE': pvalf,'ENRICHED_GENES':','.join(xl),'TERM_GENE_SET_SIZE':m})

    term_data_df = pd.DataFrame(term_data)
    if len(term_data_df) > 0:
        corrected_p = multipletests(term_data_df.PVALUE.tolist(), method='fdr_bh', is_sorted=False, returnsorted=False)
        term_data_df['QVALUE'] = list(corrected_p[1])
        return term_data_df[term_data_df['QVALUE'] < alpha]
    else:
        return pd.DataFrame()
        
def generate_gitools_matrix(cluster_genes,fnam):
    allclustergenes = []
    cluster = []
    for k,v in cluster_genes.items():
        allclustergenes += v
        cluster.append(k)
    allclustergenes = list(set(allclustergenes))
    fout = open('../results/'+fnam.split('.')[0]+'.bdm','w')
    fout.writelines('GENE\t'+'\t'.join(cluster)+'\n')
    for gene in allclustergenes:
        fout.writelines(str(gene)+'\t')
        vs = []
        for c in cluster:
            if gene in cluster_genes[c]:
                vs.append('1')
            else:
                vs.append('0')
        fout.writelines(('\t').join(vs)+'\n')
    fout.close()


def enrichment_all_groups(genesets,xmin,filt,allgenes,cluster_genes,Pvalue_thr):
    allenrichment_results = pd.DataFrame()
    for cluster, cluster_genes_l in cluster_genes.items():
        print('Working on....',cluster)
        enrichment_results = functional_enrichment(term_genes=genesets, N=len(allgenes), k_l=cluster_genes_l, alpha=Pvalue_thr, xmin=xmin)
        enrichment_results['GROUP'] = cluster
        allenrichment_results = allenrichment_results.append(enrichment_results)

    if filt != 'No_file':
        gos_sel = pd.read_excel('../data/'+filt,header=None,names=['TERM'])['TERM'].tolist()
        return allenrichment_results[allenrichment_results['TERM'].isin(gos_sel)]
    else:
        return allenrichment_results

def genes_enriched_per_group(allenrichment_results,out_ex):	
    try:
        os.stat('../results/'+out_ex)
        os.rmdir('../results/'+out_ex)
        os.mkdir('../results/'+out_ex)
    except:
        os.mkdir('../results/'+out_ex) 

    for term,termdf in allenrichment_results.groupby('TERM'):
        allgenes = list(set([gg for g in termdf.ENRICHED_GENES.tolist() for gg in g.split(',')]))
        genes_l = []
        for gene in allgenes:
            d = {} 
            for cluster, clusterdf in termdf.groupby('GROUP'):
                if gene in clusterdf.ENRICHED_GENES.tolist()[0].split(','):
                    val = 1
                else:
                    val = 0
                d[cluster] = val
            d['GENE'] = gene
            genes_l.append(d)
        term_genes = pd.DataFrame(genes_l)
        term_genes.to_excel('../results/'+out_ex+'/'+term+'.xls',index=False)
    return

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
##################################################################################################################

####################################### GSEA functions ###########################################################
def compute_median(row):
    row['rest'] = np.median([row[e] for e in row])
    return row

#Paired combinations
def paired_combinations(exp,db_name,analysis_name,permutations_number,output_plots):
    combi_cols = combinations([column for column in exp.columns.tolist() if column!='GENE'],2)
    for pair in combi_cols:
        print('---> working on combination:',' vs '.join(pair))
        exp_c = exp[['GENE',pair[0],pair[1]]].copy(deep=True)
        exp_c[pair[0]+'_log2'] = exp_c[pair[0]].apply(lambda x:compute_log2(x))
        exp_c[pair[1]+'_log2'] = exp_c[pair[1]].apply(lambda x:compute_log2(x))
        print('\t',len(exp_c[(exp_c[pair[0]+'_log2']==-666)|(exp_c[pair[1]+'_log2']==-666)]),' genes discarded for having null (0 or <0) values in one or both conditions')
        exp_c = exp_c[(exp_c[pair[0]+'_log2']!=-666)&(exp_c[pair[1]+'_log2']!=-666)]
        exp_c['FC'] = exp_c[pair[0]+'_log2'] - exp_c[pair[1]+'_log2']
        exp_c['RANK'] = rankdata(exp_c['FC'].tolist())
        print('\trunning GSEA for',pair,' this might take hours')
        run_GSEA_function(exp_c[['GENE','RANK']],db_name,analysis_name+'_'+pair[0]+'_vs_'+pair[1],permutations_number,output_plots) #run GSEA
    return

#One againts all the others
def one_againts_all_combinations(exp,db_name,analysis_name,permutations_number,output_plots):
    for c in exp.columns.tolist():
        if c!= 'GENE':
            print('---> working on combination:',c,'vs the rest')
            exp = exp.dropna()
            exp_c = exp[['GENE',c]].copy(deep=True)
            medians = [np.median([row[c2] for c2 in exp.columns.tolist() if c2 not in [c,'GENE']]) for row in exp[[c2 for c2 in exp.columns.tolist() if c2 not in [c,'GENE']]].to_dict(orient='records')]
            exp_c['rest'] = medians
            exp_c[c+'_log2'] = exp_c[c].apply(lambda x:compute_log2(x))
            exp_c['rest_log2'] = exp_c['rest'].apply(lambda x:compute_log2(x))
            print('\t',len(exp_c[(exp_c[c+'_log2']==-666)|(exp_c['rest_log2']==-666)]),' genes discarded for having null (0 or <0) values in one or both conditions')
            exp_c = exp_c[(exp_c[c+'_log2']!=-666)&(exp_c['rest_log2']!=-666)]
            exp_c['FC'] = exp_c[c+'_log2'] - exp_c['rest_log2']
            exp_c['RANK'] = rankdata(exp_c['FC'].tolist())
            exp_c[['GENE',c,'rest',c+'_log2','rest_log2','FC']].to_excel('../results/'+'FCvalues_'+analysis_name+'_'+c+'_vs_the_rest.xls',index=False)
            print('\trunning GSEA for',c,' vs the rest, this might take hours')
            run_GSEA_function(exp_c[['GENE','RANK']],db_name,analysis_name+'_'+c+'_vs_the_rest',permutations_number,output_plots) #run GSEA
    return

#GSEA function
def run_GSEA_function(genes_rankval,db_name,analysis_name,permutations_number,output_plots):
    #Run GSEA
    pre_res = gp.prerank(rnk=genes_rankval[['GENE','RANK']],
                             gene_sets=db_name,
                             outdir='../results/'+analysis_name,format='png',
                             graph_num=output_plots,
                             permutation_num=permutations_number,
                             weighted_score_type=0)
    results = pd.DataFrame(pre_res.res2d)
    results = results.reset_index()
    #Dump a summary of the results into an Excel file
    results.to_excel('../results/'+analysis_name+'_GSEAresults_SUMMARY.xls',index=False)

def compute_log2(x):
    if x <= 0:
        return -666
    else:
        return np.log2(float(x))



