import numpy as np
import pandas as pd
import ngsfragments as ngs
import os


# Local imports
from ...core.core import cfDNA


def write_seg(cfdna_object: cfDNA,
              prefix: str) -> None:
    """
    Write seg file

    Parameters
    ----------
        cfdna_object : cfDNA
            cfDNA object
        prefix : str
            Prefix for output file

    Returns
    -------
        None

    """
    
    # Write seg file
    seg_fn = prefix + ".seg"
    sample = os.path.split(prefix)[-1]
    ngs.segment.cnv_utilities.write_seg_file(cfdna_object,
                                             seg_fn,
                                             sample)
    
    return None


def write_anno(cfdna_object: cfDNA,
               prefix: str) -> None:
    """
    Write annotation file

    Parameters
    ----------
        cfdna_object : cfDNA
            cfDNA object
        prefix : str
            Prefix for output file

    Returns
    -------
        None

    
    """

    # Write annotation file
    cfdna_object.anno.engine.df.to_csv(prefix+"_metrics.txt", header=True, index=True, sep="\t")

    return None


def write_anno_seg(cfdna_object: cfDNA,
                   prefix: str,
                   bin_size: int) -> None:
    """
    Write annotationed seg file

    Parameters
    ----------
        cfdna_object : cfDNA
            cfDNA object
        prefix : str
            Prefix for output file
        bin_size : int
            Bin size
    
    Returns
    -------
        None
    """

    # Write annotation file
    sample = os.path.split(prefix)[-1]
    df = cfdna_object.obs_intervals[sample]["cnv_segments"].df
    df.loc[:,"chrom"] = cfdna_object.obs_intervals[sample]["cnv_segments"].index.labels
    df.loc[:,"start"] = cfdna_object.obs_intervals[sample]["cnv_segments"].index.starts
    df.loc[:,"end"] = cfdna_object.obs_intervals[sample]["cnv_segments"].index.ends
    df.loc[:,"sample"] = sample
    df.loc[:,"n_bins"] = ((df.loc[:,"end"].values - df.loc[:,"start"].values) / bin_size).astype(int)
    df = df.loc[:,['sample', 'chrom', 'start', 'end', 'copy_number', 'event', 'subclone_status',
                    'logR_Copy_Number', 'Corrected_Copy_Number', 'Corrected_Call', 'var', 'n_bins', 'median']]
    df.to_csv(prefix+"_seg_annotations.seg", header=True, index=False, sep="\t")

    # Drop columns
    df.drop(columns=["chrom", "start", "end"], inplace=True)

    return None