from typing import Dict, List, Tuple, Any
import pandas as pd
import numpy as np
import os
from intervalframe import IntervalFrame
from ngsfragments import Fragments
from ailist import LabeledIntervalArray

# Local imports
from ...core.core import cfDNA
from .wps import wps, coverage


def gene_sum(scores: Dict[str, pd.Series],
             gene_windows: IntervalFrame) -> pd.Series:
    """
    Sum the score values in a window

    Parameters
    ----------
        scores : Dict[str, pd.Series]
            Scores for each chromosome
        gene_windows : IntervalFrame
            Gene associated intervals

    Returns
    -------
        total_sums : pd.Series
            Total sums
    """

    # Initialize results
    chroms = list(scores.keys())
    genes = pd.unique(gene_windows.loc[chroms,"Gene"].values)
    total_sums = pd.Series(np.zeros(len(genes)),
                           index=genes)

    # Iterate over chromosomes
    for chrom in scores:
        chrom_intervals = gene_windows.loc[chrom,:]
        sum_scores = pd.Series(np.zeros(chrom_intervals.shape[0]),
                               index=chrom_intervals.df.loc[:,"Gene"].values)
        
        # Iterate over intervals
        for i, interval in enumerate(chrom_intervals.index):
            sum_scores.values[i] = np.sum(scores[chrom].loc[interval.start:interval.end].values)

        # Sum sums
        summarized_scores = sum_scores.groupby(by=chrom_intervals.df.loc[:,"Gene"].values).sum()
        total_sums.loc[summarized_scores.index] = summarized_scores.values

    return total_sums


def score_gene_sum(frags: Fragments,
                   gene_windows: IntervalFrame,
                   score_method: str = "wps",
                   adjust_lengths: bool = True,
                   protection: int = 120,
                   min_length: int = 120,
                   max_length: int = 210,
                   normalize: bool = True,
                   norm_method: str = "mean",
                   min_gene_length: int = 0,
                   max_gene_length: int = 100000) -> pd.Series:
    """
    Sum scores in all gene_windows
    """

    # Initialize results
    genes = pd.unique(gene_windows.df.loc[:,"Gene"].values)
    results = pd.Series(np.zeros(len(genes)),
                        index=genes)

    # Iterate over chroms
    chroms = set(frags.frags.unique_labels) & set(gene_windows.index.unique_labels)
    for chrom in chroms:

        # Calculate scores
        if score_method == "wps":
            scores = wps(frags, chrom, protection, min_length, max_length, normalize, norm_method)
        elif score_method == "cov":
            scores = frags.coverage(chrom, min_length=min_length, max_length=max_length)
        else:
            raise ModuleNotFoundError("Provided method is not valid [wps, cov]")

        # Sum scores
        sums = gene_sum(scores, gene_windows)
        results.loc[sums.index] = sums.values

    if adjust_lengths:
        results = length_adjust(results, gene_windows, min_length=min_gene_length, max_length=max_gene_length)

    return results


def exon_count(intervals: LabeledIntervalArray,
               gene_windows: IntervalFrame,
               adjust_lengths: bool = True) -> pd.Series:
    """
    Count n hits in a window
    """

    # Count number of hits in gene_windows
    nhits = intervals.nhits_from_LabeledIntervalArray(gene_windows.index)
    summarized_nhits = pd.Series(nhits).groupby(by=gene_windows.df.loc[:,"Gene"].values).sum()

    # Adjust lengths
    if adjust_lengths:
        summarized_nhits = length_adjust(summarized_nhits, gene_windows)

    return summarized_nhits


def length_adjust(values: pd.Series,
                  gene_windows: IntervalFrame,
                  min_length: int = 0,
                  max_length: int = 1000000) -> pd.Series:
    """
    Divide by gene length

    Parameters
    ----------
        values: pd.Series
            Scores
        gene_windows: IntervalFrame
            Gene windows
        min_length: int = 0
            Minimum gene length to include
        max_length: int = 100000
            Maximum gene length to include

    Returns
    -------
        adjusted : pd.Series
            Adjusted scores
    """

    # Find lengths
    starts = gene_windows.index.extract_starts()
    ends = gene_windows.index.extract_ends()
    l = pd.Series(ends-starts).groupby(by=gene_windows.df.loc[:,"Gene"].values).sum()

    # Adjust by dividing
    adjusted = values / l.loc[values.columns.values]

    # Filter by length
    chosen = np.logical_and(l.loc[values.columns.values] > min_length,
                            l.loc[values.columns.values] < max_length)
    adjusted = adjusted.loc[:,chosen.values]

    return adjusted


def call_wps_sum(cfdna_object: cfDNA,
                 frags: Fragments,
                 windows: List[str] = ["exons", "genes", "tss"],
                 genome: str = "hg19",
                 protection: int = 60,
                 min_length: int = 100,
                 max_length: int = 210,
                 normalize: bool = True,
                 norm_method: str = "mean",
                 adjust_lengths: bool = True,
                 use_coverage: bool = False,
                 params: Dict[str,Any] | None = None) -> None:
    """
    """

    # Determine name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Log frags
    cfdna_object.log_fragments(frags)

    # Default parameters
    default_params = {"exons": {"upstream": 5000},
                      "genes": {"upstream": 5000,
                                "gene_type": "protein_coding"},
                      "tss": {"upstream": 1000,
                              "downstream": 1000}}
    if params is None:
        params = default_params

    # Get windows
    if genome == "hg19":
        from hg19genome import get_genes as GenomeCaller
    else:
        raise NotImplementedError("Only hg19 is supported currently")
    gene_windows = {w:GenomeCaller(w, **params[w]) for w in windows}

    # Initialize results
    results = {}
    for w in windows:
        genes = pd.unique(gene_windows[w].df.loc[:,"Gene"].values)
        genes.sort()
        results[w] = pd.DataFrame(np.zeros((1, len(genes))),
                        columns=genes, index=[file_name])


        # Iterate over chroms
    for chrom in frags.chroms:
        if use_coverage:
            scores = coverage(frags, chrom, min_length, max_length, normalize, norm_method)
        else:
            scores = wps(frags, chrom, protection, min_length, max_length, normalize, norm_method)
        # Iterate over windows
        for w in windows:
            # Sum scores
            sums = gene_sum(scores, gene_windows[w])
            results[w].loc[file_name, sums.index] = sums.values

    if adjust_lengths:
        for w in windows:
            results[w] = length_adjust(results[w], gene_windows[w])

    for w in windows:
        cfdna_object.add_obs_values(file_name, w+"_sum_wps", results[w])

    return None