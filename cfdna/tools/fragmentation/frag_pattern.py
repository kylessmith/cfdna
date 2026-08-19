import numpy as np
import pandas as pd
import ngsfragments as ngs
import os
from intervalframe import IntervalFrame

# Local imports
from ...core.core import cfDNA


def fragment_profile(frags: ngs.Fragments,
                     bin_size: int = 100000,
                     genome_version: str = "hg38",
                     bias_correct: bool = True,
                     verbose: bool = False):
   """
   Determine fragment profile (ratio of large fragments to small fragments)

   Parameters
   ----------
      frags : ngsfragments.fragments
      chrom
            str
      bins_size
            int
      genome
            str
      bin_bias_h5_fn
            str
      bias_correct
            bool
      smooth
            bool

   Returns
   -------
      fragment_profile
            pandas.DataFrame
   """

   # Calculate coverage of small and large fragments
   small_bin_coverage = ngs.coverage.bin_counts(frags, bin_size=bin_size, min_length=60, max_length=160)
   large_bin_coverage = ngs.coverage.bin_counts(frags, bin_size=bin_size, min_length=161, max_length=250)

   # Correct coverage
   if bias_correct:
      cnv = ngs.segment.CNVcaller(cnv_binsize=bin_size, genome_version=genome_version)
      small_bin_coverage = cnv.correct_bins(small_bin_coverage)
      large_bin_coverage = cnv.correct_bins(large_bin_coverage)

   # Calculate fragment profile
   fragment_profile = small_bin_coverage.copy()
   fragment_profile.df.loc[:,"counts"] = small_bin_coverage.df.loc[:,"counts"].values - large_bin_coverage.df.loc[:,"counts"].values
   fragment_profile.df.columns = ["fragment_profile"]
   
   return fragment_profile


def fragment_entropy(frags: ngs.Fragments,
                     bin_size: int = 100000,
                     genome_version: str = "hg38",
                     min_length = 60,
                     max_length = 250,
                     n_bins = 5,
                     normal_bins: IntervalFrame = None,
                     rpca_normal: bool = False,
                     normalize_entropy: bool = True,
                     verbose: bool = False):
      """
      """
      from scipy.stats import entropy

      # Create bins
      cnv = ngs.segment.CNVcaller(cnv_binsize=bin_size, genome_version=genome_version)
      bins = cnv.bin_data(frags)
      cbins = cnv.correct_bins(bins)
      shannon_entropy = np.zeros(cbins.shape[0])
      sample = bins.df.columns.values[0]

      # Iterate over chromosomes
      for i, interval in enumerate(cbins.index):
            if verbose and i % 100 == 0: print(f"{i}/{cbins.shape[0]}", flush=True)
            bin_frags = frags.frags.intersect(interval.start, interval.end, interval.label).filter(min_length, max_length)
            length_dist = np.zeros(max_length-min_length)
            ld = bin_frags.length_dist()[min_length:]
            length_dist[:len(ld)] = ld
            hist_lengths = length_dist.reshape(-1,n_bins).sum(axis=1)
            hist_lengths = hist_lengths / hist_lengths.sum()
            shannon_entropy[i] = entropy(hist_lengths, base=2)

      # Create DataFrame
      cbins.df.columns = ["cnv_ratio"]
      cbins.df.loc[:,"entropy"] = shannon_entropy

      # Normalize entropy
      cbins.df = cbins.df.drop(columns=["cnv_ratio"])
      if normalize_entropy:
            cbins.df.columns = ["counts"]
            cnv = ngs.segment.CNVcaller(cnv_binsize=bin_size, genome_version=genome_version, rpca_normal=rpca_normal)
            cbins = cnv.correct_bins(cbins, normal_bins=normal_bins)
            cbins.df.columns = ["entropy"]

      return cbins

