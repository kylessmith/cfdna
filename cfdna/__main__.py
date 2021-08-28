import argparse
from .commandline.commands import cfDNA_summarize, cfDNA_merge


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    # Summarize
    parser_summary = subparsers.add_parser("summarize", help="Calculate, plot, and save summary metrics")
    parser_summary.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_summary.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_summary.add_argument("--bin_size", type=int, help="Bin size to use (default=100000)", default=100000)
    parser_summary.add_argument("--genome", help="Version of genome to use (default=hg19)", default="hg19")
    parser_summary.add_argument("--proportion", type=float, help="Proportion of fragments to use (default: 1.0)", default=1.0)
    parser_summary.add_argument("--n_frags", type=int, help="Estimate of number of fragments to use (default: ALL)", default=0)
    parser_summary.add_argument("--n_jobs", type=int, help="Number of CPUs to use (default: 1)", default=1)
    parser_summary.add_argument("--min_length", type=int, help="Minimum length to consider (default: 1)", default=1)
    parser_summary.add_argument("--max_length", type=int, help="Maximum length to consider (default: 1000)", default=1000)
    parser_summary.add_argument("--segs", help="Whether to write seg files (default: False)", default=False, action="store_true")
    parser_summary.set_defaults(func=cfDNA_summarize)

    # Merge summaries
    parser_merge = subparsers.add_parser("merge", help="Calculate, plot, and save summary metrics")
    parser_merge.add_argument("--inputs", help="Directory for h5 summary files", required=True, nargs='*')
    parser_merge.add_argument("--output", help="h5 file name for output summary", required=True)
    parser_merge.set_defaults(func=cfDNA_merge)

    args = parser.parse_args()

    args.func(args)