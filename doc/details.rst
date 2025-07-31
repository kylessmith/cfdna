cfDNA-CNVcaller detailed description
====================================

Description
-----------

The cfdna package is a python package that utilizes genome-wide coverage from aligned short read sequencing files (bam) to predict tumor fraction based on the presence and amplitude of copy-number variations (CNVs) within a sample. Adjacent DNA fragment sums for non-overlapping bins are clustered into segments based

Introduction
------------

Copy number variations (CNVs) represent a significant class of genomic variations where segments of DNA are either amplified or deleted, resulting in abnormal copy numbers in these regions. This phenomenon resulting from chromosomal instability is observed across different type of tumors. In general, copy number calling involves detecting these alterations from utilizing next-generation sequencing (NGS). Such methods help identify regions of the genome that have been duplicated or lost, providing valuable insights into cancer genetics.

Liquid biopsies represent a revolutionary advancement in this field. Instead of invasive tissue sampling, liquid biopsies analyze circulating tumor DNA (ctDNA) found in bodily fluids—primarily blood. This approach enables non-invasive detection of CNAs, offering significant advantages for monitoring disease progression, treatment response, and detecting minimal residual disease.

The technical challenges in copy number calling from liquid biopsies are substantial, including low ctDNA concentrations, high background noise from non-tumor DNA, and the need for sensitive analytical methods. Despite these challenges, recent technological advances in sequencing technologies and computational algorithms have improved the accuracy and reliability of CNA detection in liquid biopsies, making them increasingly valuable in clinical oncology.

In the next section, we provide a formal description of the cfdna-CNVcaller. At a high level, the package transforms aligned read coverages into a tumor fraction prediction.

Overview of the cfdna-CNVcaller architecture
--------------------------------------------

To infer disease status from CSF-derived whole-genome sequencing data, we developed a computational framework to detect large-scale CNVs (cfdna v#). For each sample, mean fragment coverage was calculated for 1Mb non-overlapping bins throughout the genome and filtered to remove common blacklisted regions. Bin sums were then corrected for GC content and mapping using LOESS regression and normalized by taking the log2 ratio between each bin and the median bin coverage (ngsfragments v#). Copy-number segments were determined from resulting bins using a Bayesian change point approach with a heuristic probability threshold of 0.3 (bcpseg #). A probabilistic Hidden Markov Model (HMM) model was implemented in Python and integrates into downstream analysis pipelines. The resulting Python package was used to predict tumor purity and ploidy based on larger 1Mb bins (hmmcnv v#). Each chromosome is binned into 1Mb segments and assigned the median log2 ratio within each window. Pearson correlation coefficients were calculated for all matched autosomal segments between primary tumor and baseline cfDNA.

The key components of the cfdna-CNVcaller are:
	1.	Preprocessing
		a.	Normalization
		b.	Blacklisted regions
		c.	GC and mappability correction
	2.	Segmentation
		a.	Bayesian changepoint segmentation
	3.	Tumor fraction prediction
		a.	HMM model training
		b.	Segment classification
		c.	HMM model prediction
		
Design Considerations and Advantages
------------------------------------

lcWGS of cfDNA presents several analytical challenges including (a) very low coverage of sequencing; (b) absence of matched normal germline DNA; and (c) low tumor content of many cfDNA samples. Therefore, we implemented a solution that accounts for these challenges using several assumptions that will help with the analysis and interpretation.
	(1) Large-scale CNA can be detected by evaluating read coverage in large, equal-sized genomic windows (or bins).
	(2) Homozygous deletions are typically at smaller scales than the large bin sizes used here and are not considered.
	(3) Clonal copy number states should be discrete integers.
	(4) This bin size is large enough to overcome any biases related to nucleosome positioning which is at the scales of 166 bp and 332 bp.
	(5) Due to the low coverage and absence of allelic information, only one subclone is assumed to be detectable.
	(6) The last assumption is a consequence of a limitation of the algorithm to reliably and explicitly distinguish large numbers of subclones, given the low coverage and low tumor content.

Detailed description of cfdna-CNVcaller architecture
----------------------------------------------------

cfdna-CNVcaller consists of 3 main parts:
	1.	Preprocessing
	2.	Segmentation
	3.	Tumor fraction prediction

Preprocessing

The preprocessing step of the cfdna-CNVcaller extracts DNA fragment bounds from an aligned reads file, sums the number overlapping genome tiles and performs data correction computations. There are 3 parts to the preprocessing:
	1.	Normalization
	2.	Removing blacklisted regions
	3.	GC and mappability correction
Correspondingly, we compute
	Non-overlapping bin sums based on the total number of non-duplicated DNA fragments present in each bin and remove bins consisting of greater than 50% blacklisted regions from the genome. Blacklisted regions of the genome were determined using a combination of the ENCODE blacklist and regions found to be recurrently altered (>2 std.) in a at least 2 of 7 previously published panel of normal profiles (site). 

Segmentation
	
	Non-overlapping bin sums from the preprocessing step are segemented via an online bayesian changepoint algorithm per chromosome.

Tumor fraction prediction

	A hidden markov model (HMM) is fitted to the non-overlapping bin sums across the genome. 