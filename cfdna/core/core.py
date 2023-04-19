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
                 ref_genome: str = "hg19"):
        """
        Initialize cfDNA object
                 
        Parameters
        ----------
            frags : :class: `ngsfragments.ngsfragments` (default: None)
                DNA fragment intervals from BAM file
            ref_genome : str
                Reference genome
            verbose : bool
                 Print process
                 
        Returns
        -------
            None

        """

        super().__init__(anno, values, obs_values, intervals,
                         obs_intervals, obs, uns, params)
        
        self.params["ref_genome"] = ref_genome


    def log_fragments(self,
                      frags: Fragments) -> None:
        """
        Log fragments

        Parameters
        ----------
            frags : :class: `ngsfragments.ngsfragments`
                DNA fragment intervals from BAM file
        
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
        #path = os.path.normpath(frags.sam_file)
        #file_name = path.split(os.sep)[-1]
        file_name = frags.sam_file

        # Check is file_name exists
        if file_name in self.obs:
            return

        else:
            # Find all file names
            self.add_obs(file_name)

            # Record
            self.values("length_dist", frags.length_dist())
            self.add_anno("n_fragments", file_name, len(frags.frags))
            self.uns["chrom_lengths"] = frags.genome
            self.uns["chroms"] = frags.chroms


        
