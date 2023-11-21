from typing import Dict, List, Tuple
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from joblib import Parallel, delayed
from ...core.core import cfDNA
from ...utilities.h5_utilities import merge_h5
from ..fragmentation.frag_pattern import fragment_profile
from ..cnv.segmentation import call_cnv_pipline


def call_gene_wps():
    pass

def call_tss_wps():
    pass

def call_frag_profile(cfdna_object, frags, bin_size=100000, bin_bias_h5_fn=None, bias_correct=True, smooth=True, verbose=False):
    fragment_profile(cfdna_object, frags, bin_size=100000, bin_bias_h5_fn=None, bias_correct=True, smooth=True, verbose=False)

def call_cnv_pipeline(cfdna_object, frags, cnv_binsize=cnv_binsize, hmm_binsize=hmm_binsize, method=method, outlier_smooth=outlier_smooth,
                         gauss_smooth=gauss_smooth, bcp_cutoff=bcp_cutoff, normal=normal, ploidy=ploidy, estimatePloidy=estimatePloidy,
                         minSegmentBins=minSegmentBins, maxCN=maxCN, verbose=verbose):
    call_cnv_pipline(cfdna_object, frags, cnv_binsize=cnv_binsize, hmm_binsize=hmm_binsize, method=method, outlier_smooth=outlier_smooth,
                         gauss_smooth=gauss_smooth, bcp_cutoff=bcp_cutoff, normal=normal, ploidy=ploidy, estimatePloidy=estimatePloidy,
                         minSegmentBins=minSegmentBins, maxCN=maxCN, verbose=verbose)


def pipeline(bam_filename: str,
             functions: List[str],
             **kwargs: int) -> str:
    """
    Pipe functions together while optimizing for reducing
    redundent calculations

    Parameters
    ----------
        bam_filename : str
            Bam file name
        functions : List[str]
            List of functions
        **kwargs : Iterable
            Key word arguements

    Returns
    -------
        h5_file : str
            Name of output h5 file
    """

    # Function calls
    function_calls = {"gene_wps": (call_gene_wps, "wps"),
                      "tss_wps": (call_tss_wps, "wps"),
                      "fragment_profile": (call_frag_profile, "frag"),
                      "cnv": (call_cnv_pipeline, "cnv")}
    
    # Determine functions
    func_typing = {"wps": [],
                   "frag": [],
                   "cnv": []}
    for func in functions:
        func_typing[function_calls[func][1]].append(function_calls[func][0])

    # Initialize cfdna object
    cfdna_object = cfDNA()

    # Read bam file
    frags = readBam(bam_filename)

    # Iterate over wps functions
    for func in func_typing["wps"]:
        analysis_func = func_typing["wps"][func]
        
        # Iterate over chromosomes
        for chrom in frags.frags.unique_labels:
            wps_scores = wps(frags, chrom=chrom)
            analysis_func(cfdna_object, wps_scores)


def multi_apply(bam_files, function, njobs, output_h5=None, merge=True, cleanup=True, **kwargs):
    """
    """

    function_calls = {"gene_wps": call_gene_wps,
                      "tss_wps": call_tss_wps,
                      "fragment_profile": call_frag_profile,
                      "cnv": call_cnv_pipeline}

    # Call function and write h5 files
    h5_filenames = Parallel(n_jobs=njobs)(delayed(function)(bam, **kwargs) for bam in bam_files)

    # Merge
    if merge:
        if output_h5 is None:
            output_h5 = function.__name__ + "_appied.h5"
        merge_h5(h5_filenames, output_h5, progress_bar=False)

    # Clean up
    if merge and cleanup:
        for filename in h5_filenames:
            os.remove(filename)

    return None