import pandas as pd
import numpy as np
import ngsfragments as ngs
import genome_info

def read_end_motifs(frags: ngs.Fragments,
                    genome: str = "hg38",
                    verbose: bool = False) -> pd.DataFrame:
    """
    """

    # Get genome info
    genome = genome_info.GenomeInfo(genome)

    # Store motifs
    motifs = []
    for b1 in ["A", "C", "G", "T"]:
        for b2 in ["A", "C", "G", "T"]:
            for b3 in ["A", "C", "G", "T"]:
                for b4 in ["A", "C", "G", "T"]:
                    motifs.append(b1 + b2 + b3 + b4)
    conversion = {"A": "T", "C": "G", "G": "C", "T": "A"}
    motifs = pd.DataFrame(0.0, index=motifs, columns=["front", "end"])

    # Get front motifs
    chromosomes = genome["autosomes"]
    if verbose: print("Getting front motifs...", flush=True)
    for chrom in chromosomes:
        if verbose: print(chrom, flush=True)
        i = genome.interval_kmers(frags.frags.get(chrom), k=4, last_n=4)
        for m in i:
            if "N" not in m:
                motifs.loc[m, "front"] += i[m]
    
    # Get end motifs
    if verbose: print("Getting end motifs...", flush=True)
    for chrom in chromosomes:
        if verbose: print(chrom, flush=True)
        i = genome.interval_kmers(frags.frags.get(chrom), k=4, first_n=4)
        for m in i:
            if "N" not in m:
                # Convert to reverse complement
                m_rev = "".join([conversion[b] for b in m[::-1]])
                motifs.loc[m_rev, "end"] += i[m]

    return motifs