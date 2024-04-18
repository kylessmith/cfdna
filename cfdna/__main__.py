import argparse
from .commandline.commands import CNV_calling, GeneActivity


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    # Call CNVs
    parser_cnv = subparsers.add_parser("callCNVs", help="Plot and write CNV segments")
    parser_cnv.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_cnv.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_cnv.add_argument("--bin_size", type=int, help="Bin size to use (default=100000)", default=100000)
    parser_cnv.add_argument("--hmm_bin_size", type=int, help="Bin size to use (default=1000000)", default=1000000)
    parser_cnv.add_argument("--genome", help="Version of genome to use (default=hg19)", default="hg19")
    parser_cnv.add_argument("--proportion", type=float, help="Proportion of fragments to use (default: 1.0)", default=1.0)
    parser_cnv.add_argument("--n_frags", type=int, help="Estimate of number of fragments to use (default: ALL)", default=0)
    parser_cnv.add_argument("--min_length", type=int, help="Minimum length to consider (default: 1)", default=1)
    parser_cnv.add_argument("--max_length", type=int, help="Maximum length to consider (default: 1000)", default=1000)
    parser_cnv.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_cnv.add_argument("--segs", help="Whether to write seg files (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_cnv.add_argument("--cache", help="Name of h5 file to append results to", default="")
    parser_cnv.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--use_normal", help="Whether to correct using normal profiles (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--add_wps", help="Whether to add a WPS plot for TSSs (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--add_sex", help="Whether to keep sex chromosomes (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--clonal", help="Whether to predict clonality (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--anno_file", help="Whether to write text file with predictioned metrics (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--anno_segs", help="Whether to annotated seg file (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_cnv.set_defaults(func=CNV_calling)

    # Call Gebe activity
    parser_gene = subparsers.add_parser("geneActivity", help="Write gene activity")
    parser_gene.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_gene.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_gene.add_argument("--genome", help="Version of genome to use (default=hg19)", default="hg19")
    parser_gene.add_argument("--feature", help="Feature to use (default=gene)", default="gene")
    parser_gene.add_argument("--min_length", type=int, help="Minimum length to consider (default: 120)", default=120)
    parser_gene.add_argument("--max_length", type=int, help="Maximum length to consider (default: 220)", default=220)
    parser_gene.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_gene.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_gene.add_argument("--cache", help="Name of h5 file to append results to", default="")
    parser_gene.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_gene.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_gene.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_gene.set_defaults(func=GeneActivity)

    args = parser.parse_args()

    args.func(args)