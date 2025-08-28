import click
import logging
import subprocess as sp
from rich.logging import RichHandler
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
from sklearn.neighbors import KernelDensity

def cross_validate_kde(X, bandwidths, n_splits=132, dates_max=None):
    X = np.asarray(X)[:, np.newaxis]
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = []
    for bw in bandwidths:
        ll_sum = 0.0
        for train_idx, test_idx in kf.split(X):
            X_train, X_test = X[train_idx], X[test_idx]
            kde = KernelDensity(bandwidth=bw)
            X_train = np.vstack([2*dates_max - X_train, X_train, -X_train])
            kde.fit(X_train)
            ll = kde.score(X_test)
            ll_sum += ll
        avg_ll = ll_sum / n_splits
        scores.append(avg_ll)
    best_idx = np.argmax(scores)
    best_bandwidth = bandwidths[best_idx]
    return best_bandwidth, scores

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    :This script is to perform leave-one-out cross-validation of KDE bandwidth
    :Usage python loocv.py bestbwloocv WGDdate.tsv
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('kswgdtsv', type=click.Path(exists=True))
@click.option('--output', '-o', default='loo_cv_Best_WGD_KDE.pdf', show_default=True, help="file name of output")
def bestbwloocv(kswgdtsv,output):
    Root_age = 210.6938
    df = pd.read_csv(kswgdtsv,header=0,index_col=None,sep='\t').dropna().drop_duplicates(subset=['WGD ID'])
    dates = np.array(df['Consensus Mean'].to_list())
    n = len(dates)
    bandwidths = np.linspace(1, 20, 200)
    best_bw, avg_lls = cross_validate_kde(dates, bandwidths, n_splits=132, dates_max=Root_age)
    d_lls = np.diff(avg_lls)
    d_bandwidths = np.diff(bandwidths)
    slope = d_lls / d_bandwidths
    slope_threshold_alpha = 0.05
    slope_threshold = slope_threshold_alpha * (max(avg_lls) - min(avg_lls))
    k = 10
    for i in range(len(slope) - k):
        if np.all(slope[i:i+k] < slope_threshold):
            flat_start_idx = i + 1
            break
    fig,axes = plt.subplots(1,2,figsize=(16, 8))
    ax2,ax = axes.flatten()
    dates = np.asarray(dates)[:, np.newaxis]
    flatten_bw = bandwidths[flat_start_idx]
    bw_toplot= np.linspace(flatten_bw,best_bw,10)
    for bw in bw_toplot:
        kde = KernelDensity(bandwidth=bw)
        kde.fit(np.vstack([2*Root_age - dates, dates, - dates]))
        x_plot = np.linspace(0, dates.max(), 1000)[:, np.newaxis]
        log_dens = kde.score_samples(x_plot)
        density = np.exp(log_dens)
        ax.plot(x_plot[:, 0], density, color='black', lw=2, ls = '-')
        ax.set_xlabel("Age (mya)")
        ax.set_ylabel("Density")
    ax2.plot(bandwidths,avg_lls,'k-',lw=2)
    ax2.axvspan(flatten_bw,best_bw,color='gray',alpha=0.5,label='Possible bandwidth range')
    ax2.set_xlabel("Bandwidth",fontsize=10)
    ax2.set_ylabel("Average log-likelihood per test point",fontsize=10)
    ax2.legend(loc='best',frameon=False,fontsize=10)
    ax2.text(-0.05, 1.0, 'a)',transform=ax2.transAxes,fontsize=10,weight='bold')
    ax.text(-0.05, 1.0, 'b)',transform=ax.transAxes,fontsize=10,weight='bold')
    ax.invert_xaxis()
    fig.tight_layout()
    fig.savefig(output)
    plt.close()

if __name__ == '__main__':
	cli()
