## Usage python bcpa.py bcpanalysis WGDdate.tsv Maintree.txt -o result.pdf
import click
import logging
from rich.logging import RichHandler
import copy
import pandas as pd
import numpy as np
from Bio import Phylo
import pymc3 as pm
import arviz as az
from scipy import stats
import matplotlib.pyplot as plt

def getmiddlehei(bins):
    mid_values = []
    for i in range(len(bins)-1):
        mid_values.append(np.mean([bins[i+1],bins[i]]))
    return np.array(mid_values)

def kde_mode(kde_x, kde_y):
    maxy_iloc = np.argmax(kde_y)
    mode = kde_x[maxy_iloc]
    return mode, max(kde_y)

def get_totalH(Hs):
    CHF = 0
    for i in Hs: CHF = CHF + i
    return CHF

def addvvline(ax,xvalue,color,lstyle,labell):
    if labell == '': ax.axvline(xvalue,color=color, ls=lstyle, lw=1)
    else: ax.axvline(xvalue,color=color, ls=lstyle, lw=1, label='{}: {:.1f}'.format(labell,xvalue))
    return ax

def Bayesianchangpoint(ax_posterior,bins1,n1):
    age1,desity1 = getmiddlehei(bins1),n1
    distance_data = np.diff(desity1)
    distance_data = np.absolute(distance_data)
    x = age1[1:]
    y = distance_data
    y = y[x<126]
    x = x[x<126]
    assert len(x) == len(y)
    with pm.Model() as model:
        tau = pm.Uniform('tau', lower=0, upper=x.max())
        mean1 = pm.HalfNormal('mean1', sigma=1)
        mean2 = pm.HalfNormal('mean2', sigma=1)
        sigma = pm.HalfNormal('sigma', sigma=1)
        y_obs = pm.Normal('y_obs',mu=pm.math.switch(tau >= x, mean1, mean2),sigma=sigma,observed=y)
        trace = pm.sample(20000, return_inferencedata=True,tune=20000, chains=4)
        summary = az.summary(trace)
        tau_samples = trace.posterior['tau'].values.flatten()
        kde_x = np.linspace(tau_samples.min()-1,tau_samples.max()+1,num=500)
        kde_y=stats.gaussian_kde(tau_samples,bw_method='scott').pdf(kde_x)
        mode, maxim = kde_mode(kde_x, kde_y)
        Hs, Bins, patches = ax_posterior.hist(tau_samples,bins=50,color='gray', alpha=0.6, rwidth=0.8, density=False)
        CHF = get_totalH(Hs)
        scaling = CHF*1
        ax_posterior.axvline(mode, color='k', ls='--')
        ax_posterior.set_title("Posterior distribution of shift point",fontsize=15)
        ax_posterior.set_xlabel("Date (mya)",fontsize=15)
        ax_posterior.set_ylabel("Frequency",fontsize=15)
    return mode

def plotwgdcumuhistwithstempaper(ax,X,X2,label="Angiosperm WGD",color='k',bw=1,lw=3,alpha=0.6,figg=None):
    lower_stem,upper_stem = min(X2)*0.9,max(X2)*1.1
    n1, bins1, patches1 = ax.hist(X, np.arange(lower_stem,upper_stem,bw), density=True, color=color,lw=lw, histtype="step", cumulative=-1,alpha=alpha)
    inset_ax = figg.add_axes([1/2, 2/5, 2/5, 3/10])
    shiftpoint = Bayesianchangpoint(inset_ax,bins1,n1)
    ax = addvvline(ax,shiftpoint,'k','--','Estimated shift point')
    ax.legend(loc=0,fontsize=15,frameon=False)
    ax.set_ylabel("Cumulative density", fontsize = 15)
    ax.set_xlabel("Date (mya)", fontsize = 15)
    ax.set_yticks(np.arange(11)/10)
    ax.set_title("Dynamics of WGD in angiosperms",fontsize=15)

def getalllineagesage(tree):
    tree_copy = copy.deepcopy(tree)
    internal_nodes = [node for node in tree_copy.get_nonterminals()]
    root_age = tree_copy.distance(tree_copy.get_terminals()[0])
    internal_nodes_age = [root_age-tree_copy.distance(clade) for clade in internal_nodes]
    return internal_nodes_age

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    This script is to conduct the Bayesian change point analysis
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('wgdtsv', type=click.Path(exists=True))
@click.argument('tree', type=click.Path(exists=True))
@click.option('--output', '-o', default='output.pdf', show_default=True, help="file name of output")
def bcpanalysis(wgdtsv,tree,output):
    df = pd.read_csv(wgdtsv,header=0,index_col=None,sep='\t')
    df = df.dropna()
    df = df.drop_duplicates(subset=['WGD ID'])
    allinternalnodeages = getalllineagesage(Phylo.read(tree,format='newick'))
    age_n = np.array(allinternalnodeages)*100
    fig,ax = plt.subplots(1,1,figsize=(8, 4))
    plotwgdcumuhistwithstempaper(ax,df["Consensus Mean"].to_list(),age_n,label="Angiosperm WGD",color='black',bw=1,figg=fig)
    fig.tight_layout()
    fig.savefig(output)
    plt.close()

if __name__ == '__main__':
        cli()
