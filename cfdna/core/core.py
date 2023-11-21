import os
from ngsfragments import Fragments
from dataclasses import field
from projectframe import ProjectFrame, Frame, MultiFrame, MultiIntervalFrame, UnstructuredLookup, ObsFrame, ObsIntervalFrame


# Local imports

class cfDNA(ProjectFrame):
    """
    Wrapper of fragments for cell-free DNA

	:class:`~cfDNA.cfDNA` stores a fragments object

    """

    def __init__(self,
                 anno: Frame = None,
                 values: MultiFrame = None,
                 obs_values: ObsFrame = None,
                 intervals: MultiFrame  = None,
                 obs_intervals: ObsIntervalFrame  = None,
                 obs: set = None,
                 uns: UnstructuredLookup = None,
                 params: UnstructuredLookup = None,
                 genome_version: str = "hg19"):
        """
        Initialize cfDNA object
                 
        Parameters
        ----------
            frags : :class: `ngsfragments.ngsfragments` (default: None)
                DNA fragment intervals from BAM file
            genome_version : str
                Reference genome version
            verbose : bool
                 Print process
                 
        Returns
        -------
            None

        """

        super().__init__(anno, values, obs_values, intervals,
                         obs_intervals, obs, uns, params)
        
        self.params["ref_genome"] = genome_version


    def log_fragments(self,
                      frags: Fragments) -> None:
        """
        Log Fragments object to ProjectFrame

        Parameters
        ----------
            frags : Fragments
                Fragments object

        Returns
        -------
            None

        Examples
        --------
        >>> import ngsfragments
        >>> import cfDNA
        >>> frags = ngsfragments.Fragments("test.bam")
        >>> cfDNA = cfDNA.cfDNA()
        >>> cfDNA.log_fragments(frags)
        """

        # Determine name
        file_name = frags.sam_file
        sample_name = os.path.split(file_name)[-1].split(".bam")[0]

        # Check is file_name exists
        #if sample_name not in pf.obs:
        self.add_obs(sample_name)

        # Record
        if "length_dist" not in self.values.keys or pd.isnull(self.values["length_dist"].loc[sample_name,:].values).all():
            length_dist = frags.length_dist()
            length_dist.name = sample_name
            self.add_values("length_dist", length_dist)

        if "n_fragments" not in self.anno.columns or pd.isnull(self.anno.loc[sample_name, "n_fragments"]):
            self.add_anno("n_fragments", sample_name, len(frags.frags))
        
        if "chrom_lengths" not in self.uns.keys:
            self.uns["chrom_lengths"] = frags.genome.chrom_sizes
            self.uns["chroms"] = frags.chroms

        if "chroms" not in self.uns.keys:
            self.uns["chroms"] = frags.chroms

        return None
