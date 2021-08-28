import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, sem
from scipy.signal import argrelmax
from ngsfragments.plot import plot_plt2
from ngsfragments.segment.correction import gaussian_smooth
import os

# Tell matplot lib to use 'editable text'
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot_table(sam_fn, genome, bin_size, purity, ploidy, clonal, n_frags, nfr,
               median_var, font_size=14, show=True, save=None, ax=None):
    """
    Plot table of summary metrics
    """

    # Create Axes object if not given
    if ax == None:
        with sns.axes_style("ticks"):
            fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)

    # Record text
    labels = ["SAM file name:", "Tumour purity:",
              "Ploidy:", "Clonality:",
              "Number of fragments:",
              "NFR:", "Median variance:",
              "Genome:", "Bin size:"]
    numbers = [sam_fn, '{:.3f}'.format(purity),
               '{:.3f}'.format(ploidy), '{:.3f}'.format(clonal),
               str(n_frags), '{:.3f}'.format(nfr),
               '{:.3f}'.format(median_var),
               genome, str(bin_size)]
    text = list(zip(labels, numbers))

    # Table formatting
    color = "lightgrey"
    blank = "white"
    cellColours = [[color, color],
                   [blank, blank],
                   [color, color],
                   [blank, blank],
                   [color, color],
                   [blank, blank],
                   [color, color],
                   [blank, blank],
                   [color, color]]

    # Plot table
    ax.axis('tight')
    ax.axis('off')
    table = ax.table(text, loc='center', cellColours=cellColours, edges="closed")
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)

    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=80)
    if show:
        plt.show()

    return ax


def plot_genome_signal(df, signal_label, chrom_shift, line=True, plot_median=True, smooth=False, title=None,
                       ylabel=None, show=True, save=None, ax=None):
    """
    Plot a signal across the genome
    """

    # Create Axes object if not given
    if ax == None:
        with sns.axes_style("ticks"):
            fig, ax = plt.subplots(figsize=(10,5), tight_layout=True)

    # Iterate over chroms
    chrom_index = {chrom:df.loc[:,"seqnames"].values==chrom for chrom in np.unique(df.loc[:,"seqnames"].values)}
    for chrom in chrom_index:
        signal = df.loc[chrom_index[chrom], :]
        # Smooth signal
        if smooth:
            values = gaussian_smooth(signal.loc[:,signal_label].values, scale=10)
        else:
            values = signal.loc[:,signal_label].values
        shift = chrom_shift[chrom]
        if line:
            ax.plot(signal.loc[:,"start"].values+shift, values, linewidth=0.5, color="black")
        else:
            ax.scatter(signal.loc[:,"start"].values+shift, values, s=1, color="black")
        ax.axvline(x=shift, color="grey", linestyle="--")

    # Plot 0 line
    if plot_median:
        ax.axhline(y=np.median(df.loc[:,signal_label]), color="darkgrey", linewidth=0.5)
    
    # Set titles
    if title != None:
        ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel("Chromosomes", fontsize=12)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=12)
    
    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=80)
    if show:
        plt.show()


def summary(metrics, show=True, save=None):
    """
    Plot summary metrics

    Params
    ------
        metrics
            dict
        show
            bool
        save
            bool
    
    Returns
    -------
        fig
            matplotlib.pyplot.figure
    """

    # Initialize matplotlib grid
    fig = plt.figure(figsize=(10, 7), constrained_layout=True)
    gs = fig.add_gridspec(3, 3)

    # Plot table
    table_ax = fig.add_subplot(gs[0, 0])
    sam_fn = os.path.split(metrics["sam"])[-1]
    cnv_hmm_model = np.argmax(metrics.hmm_loglik.loc[:,"loglik"].values)
    purity = 1 - float(metrics.hmm_loglik.loc[cnv_hmm_model,"n_est"])
    ploidy = float(metrics.hmm_loglik.loc[cnv_hmm_model,"phi_est"])
    clonal = float(metrics.hmm_loglik.loc[cnv_hmm_model,"Frac_genome_subclonal"])

    table_ax = plot_table(sam_fn, metrics["genome"], metrics["bin_size"],
                          purity, ploidy, clonal, metrics["n_frags"], metrics["nfr"][0],
                          metrics["median_var"], font_size = 5, show=False, ax=table_ax)

    # Plot fragment length distribution
    with sns.axes_style("ticks"):
        len_ax = fig.add_subplot(gs[0, 1])
    len_ax = plot_plt.fragment_distribution(metrics["length_dist"], title=None, show=False, ax=len_ax)
    
    # Plot TSS profile
    with sns.axes_style("ticks"):
        tss_ax = fig.add_subplot(gs[0, 2])
    tss_ax = plot_plt.mean_window(metrics["tss_profile"], center_x=1000, title=None,
                                  ylabel="Normalized WPS", xlabel="Distance to TSS", show=False, ax=tss_ax)

    # Plot CNV plot
    with sns.axes_style("ticks"):
        cnv_ax = fig.add_subplot(gs[1:3,:])
    cnv_ax = plot_plt.plot_cnv(metrics, title=None, show=False, ax=cnv_ax)

    # Plot WPS peak distances
    #with sns.axes_style("ticks"):
        #peak_dist_ax = fig.add_subplot(gs[2,:])
    #peak_dist_ax = plot_genome_signal(metrics["peak_dist"], "peak_dist", metrics.chrom_shift, line=True, smooth=True, ax=peak_dist_ax, ylabel="WPS peak distances", show=False)

    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=80)
    if show:
        plt.show()

    return fig