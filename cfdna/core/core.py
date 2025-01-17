import os
from ngsfragments import Fragments
from dataclasses import field
from projectframe import ProjectFrame, Frame, MultiFrame, MultiIntervalFrame, UnstructuredLookup, ObsFrame, ObsIntervalFrame


# Local imports

class cfDNA(ProjectFrame):
    """
    Wrapper of ProjectFrame for cell-free DNA

	:class:`~cfdna.cfDNA` stores a fragments object

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
                 genome_version: str = "hg19",
                 pf: ProjectFrame = None):
        """
        Initialize cfDNA object
                 
        Parameters
        ----------
            anno : Frame
                Annotation table
            values : MultiFrame
                Data table
            obs_values : ObsFrame
                Observed data table
            intervals : MultiFrame
                Interval table
            obs_intervals : ObsIntervalFrame
                Observed interval table
            obs : set
                Set of observed samples
            uns : UnstructuredLookup
                Unstructured data
            params : UnstructuredLookup
                Parameters
            genome_version : str
                Reference genome version (default: "hg19")
            pf : ProjectFrame
                Inpute projectframe
                 
        Returns
        -------
            None

        """

        if pf is None:
            super().__init__(anno, values, obs_values, intervals,
                            obs_intervals, obs, uns, params)
        else:
            super().__init__(pf.anno, pf.values, pf.obs_values, pf.intervals,
                            pf.obs_intervals, pf.obs, pf.uns, pf.params)
        
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
        >>> import ngsfragments as ngs
        >>> import cfdna
        >>> frags = ngs.io.read_sam("test.bam", nthreads=3, genome_version="hg19")
        >>> cfDNA = cfdna.cfDNA()
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
