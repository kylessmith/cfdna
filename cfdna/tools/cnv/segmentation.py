import ngsfragments as ngs
import os
import numpy as np
from typing import List
from ngsfragments import Fragments
from projectframe import ProjectFrame
import ngsfragments as ngs

# Local imports
from ...core.core import cfDNA


def call_cnvs(cfdna_object: cfDNA,
                data: List[str] | str | Fragments,
                cnv_binsize: int = 100000,
                hmm_binsize: int = 1000000,
                genome_version: str = "hg19",
                nthreads: int = 1,
                method: str = "online_both",
                bcp_cutoff: float = 0.3,
                hazard: int = 100,
                shuffles: int = 5000,
                p: float = 0.00005,
                normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                ploidy: List[int] = [2],
                estimatePloidy: bool = False,
                minSegmentBins: int = 25,
                maxCN: int = 5,
                wgbs: bool = False,
                verbose: bool = False,
                **kwargs) -> ProjectFrame:
    """
    Call CNVs from fragments

    Parameters
    ----------
        cfdna_object : cfDNA
            cfDNA object
        data : str | Fragments
            Fragments object or SAM file
        cnv_binsize : int
            Bin size for CNV calling
        hmm_binsize : int
            Bin size for HMM
        genome_version : str
            Genome version
        nthreads : int
            Number of threads
        method : str
            Method for CNV calling
        bcp_cutoff : float
            BCP cutoff
        hazard : int
            Hazard
        shuffles : int
            Number of shuffles
        p : float
            P-value
        normal : list
            Normal values
        ploidy : list
            Ploidy values
        estimatePloidy : bool
            Estimate ploidy
        minSegmentBins : int
            Minimum segment bins
        maxCN : int
            Maximum copy number
        wgbs : bool
            Whole genome bisulfite sequencing
        verbose : bool
            Verbose
        **kwargs : dict
            Additional arguments
    
    Returns
    -------
        pf : ProjectFrame
            ProjectFrame
    """
    
    # Call CNVs
    cfdna_object = ngs.segment.cnv_pipeline.call_cnv_pipeline(cfdna_object,
                                                        frags,
                                                        cnv_binsize = cnv_binsize,
                                                        hmm_binsize = hmm_binsize,
                                                        genome_version = genome_version,
                                                        nthreads = nthreads,
                                                        method = method,
                                                        bcp_cutoff = bcp_cutoff,
                                                        hazard = hazard,
                                                        shuffles = shuffles,
                                                        p = p,
                                                        normal = normal,
                                                        ploidy = ploidy,
                                                        estimatePloidy = estimatePloidy,
                                                        minSegmentBins = minSegmentBins,
                                                        maxCN = maxCN,
                                                        verbose = verbose)
    
    return cfdna_object