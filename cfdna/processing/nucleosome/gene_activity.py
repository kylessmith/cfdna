from ngsfragments import Fragments
from intervalframe import IntervalFrame
from math import e
import numpy as np
import pandas as pd
import glob

# Local imports
from ...core import cfDNA
from ...io.read.read_bam import readBam


def score_gene_body(bins: IntervalFrame,
                    start: int,
                    end: int,
                    chrom: str):
    """
    """

    length = end - start
    nhits = bins.intersect(start, end, chrom).df.loc[:,"counts"].values
    #weight = (1 * (1/length)) + e**-1
    #weight = 1 + e**-1
    weight = (1 + (e**-1)) * (1/length)
    score = np.sum(nhits * weight)

    return score


def score_gene_upstream(bins: IntervalFrame,
                        start: int,
                        end: int,
                        chrom: str,
                        upstream: int = 100000):
    """
    """

    length = end - start
    end = start
    start = max(start - upstream, 0)
    nhits = bins.intersect(start, end, chrom)
    
    if nhits.shape[0] > 0:
        distance = end - nhits.ends()
        #weight = ((e**(-abs(distance) / 5000)) * (1/length)) + (e**-1)
        #weight = (e**(-abs(distance) / 5000)) + (e**-1)
        weight = ((e**(-abs(distance) / 5000)) + (e**-1)) * (1/length)
        score = np.sum(nhits.df.loc[:,"counts"].values * weight)
    else:
        score = 0

    return score


def score_gene_downstream(bins: IntervalFrame,
                          start: int,
                          end: int,
                          chrom: str,
                          downstream: int = 100000):
    """
    """

    length = end - start
    start = end
    end = end + downstream
    nhits = bins.intersect(start, end, chrom)
    
    if nhits.shape[0] > 0:
        distance = end - nhits.starts()
        #weight = ((e**(-abs(distance) / 5000)) * (1/length)) + (e**-1)
        #weight = (e**(-abs(distance) / 5000)) + (e**-1)
        weight = ((e**(-abs(distance) / 5000)) + (e**-1)) * (1/length)
        score = np.sum(nhits.df.loc[:,"counts"].values * weight)
    else:
        score = 0

    return score


def gene_activity(frags: Fragments,
                  genome: str = "hg19",
                  min_length: int = 120,
                  max_length: int = 220):
    """
    """
    # Get windows
    if genome == "hg19":
        import hg19genome as Genome
    else:
        raise NotImplementedError("Only hg19 is supported currently")

    # Get genes
    genes = Genome.get_gene_body(upstream=5000, gene_type="protein_coding")

    # Determine chromosomes
    chromosomes = frags.chroms

    # Iterate over chromosomes
    scores = pd.Series(np.zeros(genes.shape[0]),
                       index = genes.loc[:,"Gene"].values)
    j = 0
    for chrom in chromosomes:
        print(chrom)
        # Get genes on chromosome
        chrom_genes = genes.loc[chrom,:]
        # Determine nhits per bin
        bins = frags.frags.get(chrom).bin_nhits(500, min_length, max_length)

        # Create interval frame
        bins_iframe = IntervalFrame(bins[0])
        bins_iframe.df.loc[:,"counts"] = bins[1]

        for i in range(chrom_genes.shape[0]):
            score = score_gene_body(bins_iframe, chrom_genes.index[i].start, chrom_genes.index[i].end, chrom)
            score += score_gene_upstream(bins_iframe, chrom_genes.index[i].start, chrom_genes.index[i].end, chrom)
            score += score_gene_downstream(bins_iframe, chrom_genes.index[i].start, chrom_genes.index[i].end, chrom)
            scores[j] = score
            j += 1

    return scores


def multi_gene_activity(directory: str,
                        genome: str = "hg19",
                        min_length: int = 121,
                        max_length: int = 375,
                        nthreads: int = 1):
    """
    """

    files = glob.glob(directory+"*.bam")

    print(files[0])
    frags = readBam(files[0], nthreads=nthreads)
    scores = gene_activity(frags, genome=genome, min_length=min_length, max_length=max_length)
    scores = scores.to_frame().T
    scores.index = [files[0]]
    for bam in files[1:]:
        print(bam)
        frags = readBam(bam, nthreads=nthreads)
        scores.loc[bam,:] = gene_activity(frags, genome=genome, min_length=min_length, max_length=max_length)

    return scores