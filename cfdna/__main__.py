import argparse
from .commandline.commands import CNV_calling, GeneActivity, GenomeWPS, GenomeCoverage, create_pon, TFwps, WPSpeaks, NucleosomeCalling, FragLengthDist, DetectEndMotifs, FragmentLengthEntropy


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    # Fragment length entropy
    parser_entropy = subparsers.add_parser("fragEntropy", help="Plot fragment length entropy")
    parser_entropy.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_entropy.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_entropy.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_entropy.add_argument("--bin_size", type=int, help="Bin size to use (default=100000)", default=100000)
    parser_entropy.add_argument("--min_length", type=int, help="Minimum length to consider (default: 60)", default=60)
    parser_entropy.add_argument("--max_length", type=int, help="Maximum length to consider (default: 250)", default=250)
    parser_entropy.add_argument("--n_bins", type=int, help="Number of bins to use (default: 5)", default=5)
    parser_entropy.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_entropy.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_entropy.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_entropy.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_entropy.add_argument("--rpca_normal", help="Whether to use rPCA for normal correction (default: False)", default=False, action="store_true")
    parser_entropy.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_entropy.set_defaults(func=FragmentLengthEntropy)

    # Fragment end motif analysis
    parser_motif = subparsers.add_parser("endMotifs", help="Write fragment end motifs")
    parser_motif.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_motif.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_motif.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_motif.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_motif.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_motif.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_motif.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_motif.add_argument("--min_length", type=int, help="Minimum length to consider (default: 60)", default=60)
    parser_motif.add_argument("--max_length", type=int, help="Maximum length to consider (default: 250)", default=250)
    parser_motif.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_motif.set_defaults(func=DetectEndMotifs)

    # Fragment length distribution
    parser_frag = subparsers.add_parser("fragDist", help="Plot fragment length distribution")
    parser_frag.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_frag.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_frag.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_frag.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_frag.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_frag.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_frag.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_frag.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_frag.set_defaults(func=FragLengthDist)

    # Nucleosome calling
    parser_nuc = subparsers.add_parser("nucleosomeCalling", help="Call nucleosome positions")
    parser_nuc.add_argument('command', choices=['dpos', 'dregion', 'dtriple'],
                        help="'dpos': call positions, 'dregion': call regions, "
                             "'dtriple': call both")
    parser_nuc.add_argument('wig_file', help='Path to input bigWig (.bw) or wiggle (.wig) file')
    parser_nuc.add_argument('--prefix', default='',
                        help='Prefix for output files')

    # Signal parameters
    parser_nuc.add_argument('-q', '--height', type=float, default=0.1,
                        help='Minimum absolute signal threshold for calling')
    parser_nuc.add_argument('--mode', choices=['abs', 'pos', 'neg'], default='abs',
                        help="'abs': bidirectional, 'pos': positive only, "
                             "'neg': negative only")
    parser_nuc.add_argument('-a', '--step', type=int, default=0,
                        help='Step/bin size in bp (0=auto: 1bp for bigWig, '
                             'detect from header for wig)')

    # Position calling parameters
    parser_nuc.add_argument('-jw', '--width', type=int, default=40,
                        help='Sliding window width for summit detection (bp)')
    parser_nuc.add_argument('-jd', '--distance', type=int, default=100,
                        help='Minimum center-to-center distance between positions')
    parser_nuc.add_argument('-e', '--edge', type=int, default=1,
                        help='Detect position edges (1=yes, 0=no)')
    parser_nuc.add_argument('-g', '--gapfill', type=int, default=0,
                        help='Fill gaps between positions (1=yes, 0=no)')
    parser_nuc.add_argument('--fuzziness', type=int, default=0,
                        help='Compute fuzziness/positioning scores (1=yes, 0=no)')

    # Region calling parameters
    parser_nuc.add_argument('-rw', '--region_width', type=int, default=40,
                        help='Minimum region width (bp)')
    parser_nuc.add_argument('-rd', '--region_dis', type=int, default=3000,
                        help='Merge distance for neighboring regions (bp)')

    # Smoothing
    parser_nuc.add_argument('-z', '--smooth_width', type=int, default=0,
                        help='Smoothing window width (bp), 0=no smoothing')

    # Differential analysis
    parser_nuc.add_argument('--wig2', default=None,
                        help='Second bigWig/wiggle file for differential analysis')

    # Output
    parser_nuc.add_argument('-o', '--out', default='result',
                        help='Output directory name')
    
    parser_nuc.set_defaults(func=NucleosomeCalling)
    

    # Call CNVs
    parser_cnv = subparsers.add_parser("callCNVs", help="Plot and write CNV segments")
    parser_cnv.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_cnv.add_argument("--pon", help="Panel of normals parquet generated from createPON", default="")
    parser_cnv.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_cnv.add_argument("--bin_size", type=int, help="Bin size to use (default=100000)", default=100000)
    parser_cnv.add_argument("--hmm_bin_size", type=int, help="Bin size to use (default=1000000)", default=1000000)
    parser_cnv.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_cnv.add_argument("--proportion", type=float, help="Proportion of fragments to use (default: 1.0)", default=1.0)
    parser_cnv.add_argument("--n_frags", type=int, help="Estimate of number of fragments to use (default: ALL)", default=0)
    parser_cnv.add_argument("--min_length", type=int, help="Minimum length to consider (default: 1)", default=1)
    parser_cnv.add_argument("--max_length", type=int, help="Maximum length to consider (default: 1000)", default=1000)
    parser_cnv.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_cnv.add_argument("--segs", help="Whether to write seg files (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_cnv.add_argument("--n_mads", type=float, help="Number of MADs to merge segments (default: 1.4826)", default=1.4826)
    parser_cnv.add_argument("--normal_mads", type=float, help="Number of MADs to call normal (default: 1.4826)", default=1.4826)
    parser_cnv.add_argument("--bcp_cutoff", type=float, help="Cutoff for bayesian changpoint segmentation (default: 0.3)", default=0.3)
    parser_cnv.add_argument("--cache", help="Name of h5 file to append results to", default="")
    parser_cnv.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--use_normal", help="Whether to correct using normal profiles (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--add_wps", help="Whether to add a WPS plot for TSSs (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--add_sex", help="Whether to keep sex chromosomes (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--clonal", help="Whether to predict clonality (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--anno_file", help="Whether to write text file with predictioned metrics (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--anno_segs", help="Whether to annotated seg file (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--bins_file", help="Whether to write bins file (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--rpca_normal", help="Whether to use rPCA for normal correction (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--nanopore", help="Whether input bam is from nanopore (default: False)", default=False, action="store_true")
    parser_cnv.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_cnv.set_defaults(func=CNV_calling)

    # Call Gene activity
    parser_gene = subparsers.add_parser("geneActivity", help="Write gene activity")
    parser_gene.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_gene.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_gene.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_gene.add_argument("--feature", help="Feature to use (default=gene)", default="gene")
    parser_gene.add_argument("--protection", type=int, help="Protection winodw size (default: 120)", default=120)
    parser_gene.add_argument("--min_length", type=int, help="Minimum length to consider (default: 120)", default=120)
    parser_gene.add_argument("--max_length", type=int, help="Maximum length to consider (default: 220)", default=220)
    parser_gene.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_gene.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_gene.add_argument("--cache", help="Name of h5 file to append results to", default="")
    parser_gene.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_gene.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_gene.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_gene.set_defaults(func=GeneActivity)

    # Call Genome WPS
    parser_wps = subparsers.add_parser("genomeWPS", help="Write genome WPS")
    parser_wps.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_wps.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_wps.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_wps.add_argument("--protection", type=int, help="Protection winodw size (default: 120)", default=120)
    parser_wps.add_argument("--min_length", type=int, help="Minimum length to consider (default: 120)", default=120)
    parser_wps.add_argument("--max_length", type=int, help="Maximum length to consider (default: 220)", default=220)
    parser_wps.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_wps.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_wps.add_argument("--center", help="Way to center", default="median")
    parser_wps.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_wps.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_wps.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_wps.set_defaults(func=GenomeWPS)

    # Call Genome coverage
    parser_cov = subparsers.add_parser("genomeCov", help="Write genome coverage")
    parser_cov.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_cov.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_cov.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_cov.add_argument("--min_length", type=int, help="Minimum length to consider (default: 1)", default=1)
    parser_cov.add_argument("--max_length", type=int, help="Maximum length to consider (default: 1000)", default=1000)
    parser_cov.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_cov.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_cov.add_argument("--center", help="Way to center", default="median")
    parser_cov.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_cov.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_cov.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_cov.set_defaults(func=GenomeCoverage)

    # Call create PON
    parser_pon = subparsers.add_parser("createPON", help="Create panel of normals")
    parser_pon.add_argument("--dir", help="Directory with BAM files", required=True)
    parser_pon.add_argument("--out", help="Ouput file name", required=True)
    parser_pon.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_pon.add_argument("--bin_size", type=int, help="Bin size (default: 10000)", default=10000)
    parser_pon.add_argument("--min_length", type=int, help="Minimum length to consider (default: 1)", default=1)
    parser_pon.add_argument("--max_length", type=int, help="Maximum length to consider (default: 1000)", default=1000)
    parser_pon.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_pon.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_pon.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_pon.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_pon.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_pon.set_defaults(func=create_pon)

    # Call TF WPS
    parser_wps = subparsers.add_parser("TFwps", help="Write TF WPS")
    parser_wps.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_wps.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_wps.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_wps.add_argument("--protection", type=int, help="Protection winodw size (default: 60)", default=60)
    parser_wps.add_argument("--min_length", type=int, help="Minimum length to consider (default: 60)", default=60)
    parser_wps.add_argument("--max_length", type=int, help="Maximum length to consider (default: 220)", default=220)
    parser_wps.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_wps.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_wps.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_wps.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_wps.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_wps.set_defaults(func=TFwps)

    # Call nucleosomes
    parser_nuc = subparsers.add_parser("WPSpeaks", help="Write WPS peaks")
    parser_nuc.add_argument("--bam", help="BAM file", required=True, nargs='*')
    parser_nuc.add_argument("--prefix", help="Prefix for ouput files", default="")
    parser_nuc.add_argument("--genome", help="Version of genome to use (default=hg38)", default="hg38")
    parser_nuc.add_argument("--protection", type=int, help="Protection winodw size (default: 120)", default=120)
    parser_nuc.add_argument("--min_length", type=int, help="Minimum length to consider (default: 120)", default=120)
    parser_nuc.add_argument("--max_length", type=int, help="Maximum length to consider (default: 220)", default=220)
    parser_nuc.add_argument("--mapq", type=int, help="Mapping quality cutoff (default: 25)", default=25)
    parser_nuc.add_argument("--nthreads", type=int, help="Number of threads to use (default=1)", default=1)
    parser_nuc.add_argument("--center", help="Way to center", default="median")
    parser_nuc.add_argument("--single", help="Whether to reads are single ended (default: False)", default=False, action="store_true")
    parser_nuc.add_argument("--qcfail", help="Whether to remove qcfail flag (default: False)", default=False, action="store_true")
    parser_nuc.add_argument("--verbose", help="Whether to be verbose (default: False)", default=False, action="store_true")
    parser_nuc.set_defaults(func=WPSpeaks)

    args = parser.parse_args()

    args.func(args)