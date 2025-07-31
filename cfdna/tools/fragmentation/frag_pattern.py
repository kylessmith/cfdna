import numpy as np
import pandas as pd
import ngsfragments as ngs
import os

# Local imports
from ...core.core import cfDNA


def fragment_profile(cfdna_object: cfDNA,
                     frags: ngs.Fragments,
                     bin_size: int = 100000,
                     bin_bias_h5_fn: str = None,
                     bias_correct: bool = True,
                     smooth: bool = True,
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

   # Determine file name
   path = os.path.normpath(frags.sam_file)
   file_name = path.split(os.sep)[-1]

   # Check if file was previously annotated
   cfdna_object.log_fragments(frags)

   # Calculate coverage of small and large fragments
   small_bin_coverage = frags.bin_counts(bin_size=bin_size, min_length=100, max_length=150)
   large_bin_coverage = frags.bin_counts(bin_size=bin_size, min_length=151, max_length=220)

   # Correct coverage
   if bias_correct:
      small_bin_coverage = ngs.segment.correction.correct(small_bin_coverage, genome=cfdna_object.ref_genome, bin_size=bin_size, verbose=verbose)
      large_bin_coverage = ngs.segment.correction.correct(large_bin_coverage, genome=cfdna_object.ref_genome, bin_size=bin_size, verbose=verbose)
      small_bin_coverage = small_bin_coverage.loc[:,["corrected_counts"]]
      small_bin_coverage.df.columns = [file_name]
      large_bin_coverage = large_bin_coverage.loc[:,["corrected_counts"]]
      large_bin_coverage.df.columns = [file_name]
   else:
      small_bin_coverage.df.columns = [file_name]
      large_bin_coverage.df.columns = [file_name]

   # Calculate fragment profile
   fragment_profile = small_bin_coverage.copy()
   fragment_profile.df.loc[:,file_name] = small_bin_coverage.df.loc[:,file_name].values / large_bin_coverage.df.loc[:,file_name].values

   # Correct regions with 0 for both
   fragment_profile.df.loc[pd.isnull(fragment_profile.loc[:,file_name].values), file_name] = 1.0
   
   # Append bins
   cfdna_object.add_intervals("frag_profile", fragment_profile)


def fragment_entropy(frags: ngs.Fragments,
                     bin_size: int = 100000,
                     min_length = 100,
                     max_length = 240,
                     n_bins = 5):
      """
      """
      from scipy.stats import entropy

      # Create bins
      cnv = ngs.segment.CNVcaller(cnv_binsize=bin_size)
      bins = cnv.bin_data(frags)
      cbins = cnv.correct_bins(bins)
      shannon_entropy = np.zeros(cbins.shape[0])
      sample = bins.df.columns.values[0]

      # Iterate over chromosomes
      for i, interval in enumerate(cbins.index):
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

      return cbins

