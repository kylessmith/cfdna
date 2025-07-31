Tutorial
=========

.. code-block:: python

	from cfdna import cfDNA, cfDNA_summary
	from cfdna.plot import plot_plt
	import numpy as np
	
	# Read sam file
	frags = cfdna.io.read_sam("test_bam.bam", ntheads = 3, genome_version="hg19")
	
	# Create cfDNA object
	cf = cfdna.cfDNA(genome_version="hg19")
	
	# Determine CNVs
	cf = cfdna.tools.call_cnvs(cf, genome_version="hg19")
	
	# Determine window protection score (WPS)
	cf = cfdna.tools.wps_windows(cf, frags, sample="test")
	
	# Write summary
	cf.to_h5("test.h5")
	
	# Read summary
	cf = cfdna.io.read_h5("test.h5")
	
	# Plot
	cfdna.plot.plot_plt.plot_cnv_summary(cf, sample="test", add_wps=True, show=True, save="test_cnvs.pdf")


Commandline
-----------

.. code-block:: bash

	$ python -m cfdna createPON --dir ./ --out hg19_pon.parquet --genome hg19 --bin_size 10000
	
Note it is alway better to make PON files with smaller bins than you plan on using because they can always be merged

.. code-block:: bash
	
	$ python -m cfdna createPON -h
	
	usage: __main__.py createPON [-h] --dir DIR --out OUT [--genome GENOME] [--bin_size BIN_SIZE]
	                             [--min_length MIN_LENGTH] [--max_length MAX_LENGTH] [--mapq MAPQ] [--nthreads NTHREADS]
	                             [--single] [--qcfail] [--verbose]

	options:
	  -h, --help            show this help message and exit
	  --dir DIR             Directory with BAM files
	  --out OUT             Ouput file name
	  --genome GENOME       Version of genome to use (default=hg38)
	  --bin_size BIN_SIZE   Bin size (default: 10000)
	  --min_length MIN_LENGTH
	                        Minimum length to consider (default: 1)
	  --max_length MAX_LENGTH
	                        Maximum length to consider (default: 1000)
	  --mapq MAPQ           Mapping quality cutoff (default: 25)
	  --nthreads NTHREADS   Number of threads to use (default=1)
	  --single              Whether to reads are single ended (default: False)
	  --qcfail              Whether to remove qcfail flag (default: False)
	  --verbose             Whether to be verbose (default: False)


.. code-block:: bash

	$ python -m cfdna callCNVs --bam test.bam --genome hg38 --nthreads 3 --segs --add_wps --anno_segs --anno_file
	
.. code-block:: bash

	$ python -m cfdna callCNVs -h
	
	usage: __main__.py callCNVs [-h] --bam [BAM ...] [--prefix PREFIX] [--bin_size BIN_SIZE]
	                            [--hmm_bin_size HMM_BIN_SIZE] [--genome GENOME] [--proportion PROPORTION]
	                            [--n_frags N_FRAGS] [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
	                            [--mapq MAPQ] [--segs] [--nthreads NTHREADS] [--cache CACHE] [--single]
	                            [--qcfail] [--use_normal] [--add_wps] [--add_sex] [--clonal] [--anno_file]
	                            [--anno_segs] [--verbose]

	options:
	  -h, --help            show this help message and exit
	  --bam [BAM ...]       BAM file
	  --pon PON             Panel of normals parquet generated from createPON
	  --prefix PREFIX       Prefix for ouput files
	  --bin_size BIN_SIZE   Bin size to use (default=100000)
	  --hmm_bin_size HMM_BIN_SIZE
	                        Bin size to use (default=1000000)
	  --genome GENOME       Version of genome to use (default=hg19)
	  --proportion PROPORTION
	                        Proportion of fragments to use (default: 1.0)
	  --n_frags N_FRAGS     Estimate of number of fragments to use (default: ALL)
	  --min_length MIN_LENGTH
	                        Minimum length to consider (default: 1)
	  --max_length MAX_LENGTH
	                        Maximum length to consider (default: 1000)
	  --mapq MAPQ           Mapping quality cutoff (default: 25)
	  --segs                Whether to write seg files (default: False)
	  --nthreads NTHREADS   Number of threads to use (default=1)
	  --cache CACHE         Name of h5 file to append results to
	  --single              Whether to reads are single ended (default: False)
	  --qcfail              Whether to remove qcfail flag (default: False)
	  --use_normal          Whether to correct using normal profiles (default: False)
	  --add_wps             Whether to add a WPS plot for TSSs (default: False)
	  --add_sex             Whether to keep sex chromosomes (default: False)
	  --clonal              Whether to predict clonality (default: False)
	  --anno_file           Whether to write text file with predictioned metrics (default: False)
	  --anno_segs           Whether to annotated seg file (default: False)
	  --bins_file           Whether to write bins file (default: False)
	  --verbose             Whether to be verbose (default: False)
	  