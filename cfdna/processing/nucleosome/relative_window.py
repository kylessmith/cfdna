import numpy as np
import pandas as pd
from ngsfragments import Fragments
from intervalframe import IntervalFrame
import glob

# Local imports
from ...core import cfDNA
from .wps import wps
from ...io.read.read_bam import readBam


def calculate_relative_nhits(frags: Fragments,
                             start: int,
                             end: int,
                             chrom: str,
                             min_length: int = 120,
                             max_length: int = 220):
    """
    """

    # Calculate gene score
    gene_score = frags.frags.nhits(start,
                                    end,
                                    chrom,
                                    min_length,
                                    max_length)
    
    # Calculate upstream score
    up_start = max(start - 10000, 0)
    up_end = start
    upstream_score = frags.frags.nhits(up_start,
                                        up_end,
                                        chrom,
                                        min_length,
                                        max_length)
    # Calculate downstream score
    down_start = end
    down_end = end + 10000
    downstream_score = frags.frags.nhits(down_start,
                                        down_end,
                                        chrom,
                                        min_length,
                                        max_length)
    
    # Calculate score
    gene_score = gene_score / (end - start)
    gene_score = max(gene_score, 0.000001)
    flank_score = (upstream_score + downstream_score) / ((up_end - up_start) + (down_end - down_start))
    flank_score = max(flank_score, 0.000001)
    score = np.log2(gene_score / flank_score)

    return score


def calculate_relative_wps(wps: pd.Series,
                             start: int,
                             end: int):
    """
    """

    # Calculate gene score
    gene_score = wps.loc[start:end].values.sum()

    # Calculate upstream score
    up_start = max(start - 10000, wps.index.values[0])
    up_end = start
    upstream_score = wps.loc[up_start:up_end].values.sum()

    # Calculate downstream score
    down_start = end
    down_end = min(end + 10000, wps.index.values[-1])
    downstream_score = wps.loc[down_start:down_end].values.sum()

    # Calculate score
    gene_score = gene_score / (end - start)
    if gene_score == 0:
        gene_score = 0.001
    flank_score = (upstream_score + downstream_score) / ((up_end - up_start) + (down_end - down_start))
    if flank_score == 0:
        flank_score = 0.001
    score = np.log2((gene_score+1) / (flank_score+1))

    return score


def gene_relative_enrichment(frags: Fragments,
                             genome: str = "hg19",
                             min_length: int = 120,
                             max_length: int = 220,
                             use_wps: bool = False,
                             protection: int = 120):
    """
    """
    # Get windows
    if genome == "hg19":
        import hg19genome as Genome
    else:
        raise NotImplementedError("Only hg19 is supported currently")

    # Get genes
    #genes = Genome.get_gene_body(upstream=5000, gene_type="protein_coding")
    genes = Genome.get_tss(upstream=1000, downstream=1000)

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
        
        # Calculate WPS
        if use_wps:
            wps_scores = wps(frags, chrom, protection, min_length, max_length)

        for i in range(chrom_genes.shape[0]):
            if use_wps:
                score = calculate_relative_wps(wps_scores[chrom],
                                                chrom_genes.index[i].start,
                                                chrom_genes.index[i].end)
            else:
                score = calculate_relative_nhits(frags,
                                                chrom_genes.index[i].start,
                                                chrom_genes.index[i].end,
                                                str(chrom),
                                                min_length,
                                                max_length)
            scores[j] = score
            j += 1

    return scores


def multi_rel_enrich(directory: str,
                        genome: str = "hg19",
                        min_length: int = 121,
                        max_length: int = 375,
                        nthreads: int = 1,
                        use_wps: bool = False):
    """
    """

    files = glob.glob(directory+"*.bam")

    print(files[0])
    frags = readBam(files[0], nthreads=nthreads)
    scores = gene_relative_enrichment(frags, genome=genome, min_length=min_length, max_length=max_length, use_wps=use_wps)
    scores = scores.to_frame().T
    scores.index = [files[0]]
    for bam in files[1:]:
        print(bam)
        frags = readBam(bam, nthreads=nthreads)
        scores.loc[bam,:] = gene_relative_enrichment(frags, genome=genome, min_length=min_length, max_length=max_length, use_wps=use_wps)

    return scores