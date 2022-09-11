#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import *
import os

import warnings
warnings.filterwarnings('ignore')


# In[2]:


SAVE = True


# In[ ]:


def draw_hists(effects, labels, name, tag, c = "#a1cf53"):
    print(name)
    sns.set(rc={"figure.figsize": (25, 5)})
    sns.set_style("ticks")
    for i, (effect, label) in enumerate(zip(effects, labels)):
        subplot(1, 3, i + 1)
        ax = sns.distplot(effect, label=label, rug=False, hist=True, color = c)
        ax.set_title(f"{tag}: {label}")
 #       plt.xticks(rotation=-90)
    subplots_adjust(wspace=0.3)
    fig = plt.gcf()
    if SAVE:
        fig.savefig(f'new_images/{name.replace(" ", "_")}.png', dpi=400)
#        fig.savefig(f'images/{name.replace(" ", "_")}.pdf', dpi=30)
#    plt.show()


def draw_corr(pvals, labels, name, tag, c = "#fdb462"):
    def draw(i, lists, labels, tag, c):
        subplot(1, 3, i)
        ax = sns.scatterplot(lists[0], lists[1], color=c)
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_title(f"{tag}: Correlation = {np.corrcoef(lists[0], lists[1])[0,1]}")

    lists = [list(l) for l in pvals]
    print(name)
    sns.set(rc={"figure.figsize": (20, 5)})
    sns.set_style("ticks")
    draw(1, [lists[1], lists[2]], [labels[1], labels[2]], tag, c)
    draw(2, [lists[1], lists[0]], [labels[1], labels[0]], tag, c)
    draw(3, [lists[2], lists[0]], [labels[2], labels[0]], tag, c)

    subplots_adjust(wspace=0.3)
    fig = plt.gcf()
    if SAVE:
        fig.savefig(f'new_images/{name.replace(" ", "_")}.png', dpi=400)
#        fig.savefig(f'images/{name.replace(" ", "_")}.pdf', dpi=30)
#    plt.show()
    
def draw_scatter(lists, labels, name, tag, c = "#fdb462"):
    sns.set(rc={"figure.figsize": (6, 6)})
    sns.set_style("ticks")
    ax = sns.scatterplot(lists[0], lists[1], color=c)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(f"{tag}: correlation = {np.corrcoef(lists[0], lists[1])[0,1]}")
    fig = plt.gcf()
    if SAVE:
        fig.savefig(f'new_images/{name.replace(" ", "_")}.png', dpi=400)
#        fig.savefig(f'images/{name.replace(" ", "_")}.pdf', dpi=30)
#    plt.show()


def filter_snps(data, col, list_rsids):
    return data[data[col].isin(list_rsids)].reset_index()


dirs = [f for f in os.listdir(".") if os.path.isdir(f) and "a_" in f]
cur_dir = dirs[0]

metal_datas = []

for cur_dir in dirs:
    print(f"{cur_dir}")
    cur_files = [f for f in os.listdir(f"{cur_dir}/data")]
    metal_data_name = [
        f for f in cur_files if f[-4:] == ".TBL" and "extended_" not in f
    ][0]
    ukb_data_name = [f for f in cur_files if "gwas.imputed_v3" in f][0]
    finn_data_name = [f for f in cur_files if "hg19lifted" in f][0]

    metal_data = pd.read_csv(f"{cur_dir}/data/{metal_data_name}", sep="\t").sort_values(
        by=["MarkerName"]
    ).reset_index()
    ukb_data = pd.read_csv(f"{cur_dir}/data/{ukb_data_name}", sep="\t").sort_values(
        by=["rsid"]
    ).reset_index()
    finn_data = pd.read_csv(f"{cur_dir}/data/{finn_data_name}", sep="\t").sort_values(
        by=["rsid"]
    ).reset_index()
    
    metal_data["MAF"] = metal_data.Freq1.apply(lambda x: min(x, 1 - x))
    
    metal_data = metal_data.rename(columns = {'Effect':'beta', 'StdErr':'se', 'P-value':'pval'})
    
    if 'beta' not in metal_data.columns:
        metal_data["se"] = 1 / np.sqrt(
            2 * metal_data.MAF * (1 - metal_data.MAF) * (metal_data.Weight + metal_data.Zscore ** 2)
        )
        metal_data['beta'] = metal_data.Zscore * metal_data.se

    effects = [
        np.array(metal_data.beta),
        np.array(ukb_data.beta),
        np.array(finn_data.beta),
    ]
    stderrs = [
        np.array(metal_data.se),
        np.array(ukb_data.se),
        np.array(finn_data.se),
    ]
    pvals = [
        -np.log10(metal_data.pval),
        -np.log10(ukb_data.pval),
        -np.log10(finn_data.pval),
    ]
     
    metal_datas.append([cur_dir, metal_data])
    labels = ["metal", "ukb", "finn"]
    draw_hists(effects, labels, f"{cur_dir} effect size", cur_dir.split('__')[0], c=sns.color_palette("light:#005BBB")[-1])
    plt.clf()
    draw_hists(stderrs, labels, f"{cur_dir} stderrs", cur_dir.split('__')[0], c=sns.color_palette("light:#005BBB")[-1])
    plt.clf()
#    draw_corr(pvals, labels, f"{cur_dir} pvals corr", cur_dir.split('__')[0], c=sns.color_palette("light:#FFD500")[-1])
#    plt.clf()
    print('\n\n\n')


# In[ ]:


labels = [md[0] for md in metal_datas]
pvals = [-np.log10(md[1].pval) for md in metal_datas]
draw_corr(pvals, labels, f"methods RHEUMA_OTHER_WIDE___M06 pvals corr", '', c=sns.color_palette("light:#FFD500")[-1])
    


# In[ ]:



