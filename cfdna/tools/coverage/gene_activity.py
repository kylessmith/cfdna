import ngsfragments as ngs
from ngsfragments import Fragments
from ..core.core import cfDNA


def gene_activity(frags: Fragments,
                  cfdna_object: cfDNA,
                  genome_version: str = "hg19",
                  feature: str = "gene",
                  min_length: int = 120,
                  max_length: int = 220,
                  verbose: bool = False) -> cfDNA:
    """
    """

    # Get gene activity
    gene_activity = ngs.metrics.gene_activity(frags,
                                            genome_version=genome_version,
                                            feature=feature,
                                            min_length=min_length,
                                            max_length=max_length,
                                            verbose=verbose)
    gene_activity = ngs.metrics.correct_gene_activity(frags,
                                                     gene_activity,
                                                     correct_cnv = False,
                                                     genome_version = genome_version,
                                                     feature = feature,
                                                     verbose = verbose)

    cfdna_object.add_values("gene_activity", gene_activity)

    return cfdna_object