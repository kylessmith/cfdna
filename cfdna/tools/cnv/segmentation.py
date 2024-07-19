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
                frags: Fragments,
                cnv_binsize: int = 100000,
                hmm_binsize: int = 1000000,
                genome_version: str = "hg19",
                nthreads: int = 1,
                normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                ploidy: List[int] = [2, 3],
                estimatePloidy: bool = False,
                minSegmentBins: int = 25,
                maxCN: int = 5,
                scStates: List[int] = None,
                use_normal: bool = False,
                add_sex: bool = False,
                verbose: bool = False) -> cfDNA:
    """
    Call CNVs from fragments

    This function corrects genome-wide bin coverage for GC, mappability, and black-listed regions.
    Segments are determined using Bayesian Changepoint Segmentation and CNVs and purity  are called 
    using an HMM model. 

    Parameters
    ----------
        cfdna_object : cfDNA
            cfDNA object
        frags : Fragments
            Fragments object
        cnv_binsize : int
            Size of bins for CNV calling (default: 100000)
        hmm_binsize : int
            Size of bins for HMM (default: 1000000)
        genome_version : str
            Reference genome version (default: "hg19")
        nthreads : int
            Number of threads to use (default: 1)
        normal : List[float]
            List of normal fractions (default: [0.1, 0.25, 0.5, 0.75, 0.9])
        ploidy : List[int]
            List of ploidies (default: [2, 3])
        estimatePloidy : bool
            Estimate ploidy (default: False)
        minSegmentBins : int
            Minimum segment bins (default: 25)
        maxCN : int
            Maximum copy number (default: 5)
        scStates : List[int]
            List of states (default: None)
        use_normal : bool
            Use normal (default: False)
        add_sex : bool
            Add sex chromosomes (default: False)
        verbose : bool
            Print process (default: False)
        
    Returns
    -------
        cfdna_object : :class: `cfdna.cfDNA`
            cfDNA object

    Examples
    --------
    >>> import ngsfragments as ngs
    >>> import cfdna
    >>> frags = ngs.io.from_sam("test.bam", genome_version="hg19", nthreads=3)
    >>> cfDNA = cfdna.cfDNA()
    >>> cfDNA = cfdna.tools.cnv.call_cnvs(cfDNA,
                                      frags,
                                      cnv_binsize=100000,
                                      hmm_binsize=1000000,
                                      genome_version="hg19",
                                      nthreads=3,
                                      normal=[0.1, 0.25, 0.5, 0.75, 0.9],
                                      ploidy=[2, 3],
                                      estimatePloidy=False,
                                      minSegmentBins=25,
                                      maxCN=5,
                                      scStates=None,
                                      use_normal=False,
                                      add_sex=True)
    >>> cfDNA

    See Also
    --------
    cfdna.cfDNA : cfDNA object

    """
    
    # Call CNVs
    cfdna_object = ngs.segment.cnv_pipeline.call_cnv_pipeline(cfdna_object,
                                                            frags,
                                                            genome_version=genome_version,
                                                            cnv_binsize=cnv_binsize,
                                                            hmm_binsize=hmm_binsize,
                                                            nthreads = nthreads,
                                                            use_normal = use_normal,
                                                            keep_sex_chroms = add_sex,
                                                            normal = normal,
                                                            ploidy = ploidy,
                                                            estimatePloidy = estimatePloidy,
                                                            scStates = scStates,
                                                            minSegmentBins = minSegmentBins,
                                                            maxCN = maxCN)
    
    return cfdna_object