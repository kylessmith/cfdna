import os
import numpy as np
import pandas as pd
import h5py
import argparse
from intervalframe import IntervalFrame
from intervalframe.write.write_h5 import write_h5_intervalframe
from intervalframe.read.read_h5 import read_h5_intervalframe
from ngsfragments import fragments, utilities
from ngsfragments.segment.cnv import call_cnvs, process_cnvs
from ngsfragments.segment import cnv_utilities
from ngsfragments.segment import correction


class cfDNA(object):
    """
    Wrapper of fragments for cell-free DNA

	:class:`~cfDNA.cfDNA` stores a fragments object

    """

    def __init__(self, frags=None, cfdna_h5=None, ref_genome="hg19", verbose=False):
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

        if cfdna_h5 is None:
            # Initialize analysis properties
            self.ref_genome = ref_genome
            self.gender = {}

            self.filenames = set()
            self.chrom_lengths = None
            self.cnv_binsize = None
            self.hmm_binsize = None
            self.profile_binsize = None
            self.hmm_states = {}
            self.chroms = None

            self._anno = None
            self._intervals = {}
            self._obs_intervals = {}
            self._obs_values = {}
        else:
            pass

        # Initialize values from fragments
        if frags is not None:
            self.log_filename(frags)


    @property
    def anno(self):
        """
        pandas.DataFrame: Observation annotations
        """

        return self._anno

    @property
    def intervals(self):
        """
        dict of IntervalFrames: Shared interval values
        """

        return self._intervals


    @property
    def obs_intervals(self):
        """
        dict of IntervalFrames: Observation specific interval values
        """

        return self._obs_intervals


    @property
    def obs_values(self):
        """
        dict of pandas.DataFrame: Shared values
        """

        return self._obs_values


    def __len__(self):
        """
        Get length from self.frags
        """

        return len(self.filenames)


    def add_anno(self, key, file_name, value):
        """
        Add value to key in .anno

        Parameters
        ----------
            key : str
                Key for .anno columns
            file_name : str
                Observation file name
            value : str, float, or int
                Values to record

        Returns
        -------
            None

        """

        if self._anno is None:
            self._anno = pd.DataFrame(value, index=[file_name], columns=[key])
        self._anno.loc[file_name, key] = value


    def add_intervals(self, key, file_name, values):
        """
        Add value to key in .anno

        Parameters
        ----------
            key : str
                Key for .anno columns
            file_name : str or list-like
                Observation file name
            values : IntervalFrame
                Values to record

        Returns
        -------
            None

        """

        try:
            self._intervals[key]
            values = values.exact_match(self._intervals[key])
            self._intervals[key] = self._intervals[key].exact_match(values)
            self._intervals[key].df.loc[:,file_name] = values.loc[:,file_name].values

        except KeyError:
            self._intervals[key] = values

        
    def add_obs_intervals(self, key, file_name, values):
        """
        Add value to key in .anno

        Parameters
        ----------
            key : str
                Key for .anno columns
            file_name : str
                Observation file name
            values : IntervalFrame
                Values to record

        Returns
        -------
            None

        """

        try:
            self._obs_intervals[key]
        except KeyError:
            self._obs_intervals[key] = {}
        self._obs_intervals[key][file_name] = values


    def add_obs_values(self, key, file_name, values):
        """
        Add value to key in .anno

        Parameters
        ----------
            key : str
                Key for .anno columns
            file_name : str
                Observation file name
            values : np.array or pandas.DataFrame
                Values to record

        Returns
        -------
            None

        """

        # Check shape
        if values.shape[0] == len(self):
            pass
        elif len(values.shape) == 2:
            if values.shape[1] == len(self):
                values = values.T
        elif len(self) == 1:
            pass
        elif len(values.shape) == 1:
            pass
        else:
            raise IndexError("Shape of values cannot be matched.")

        try:
            self._obs_values[key]
            if isinstance(values, np.ndarray):
                self._obs_values[key].loc[file_name,:] = values
            elif isinstance(values, pd.DataFrame):
                self._obs_values[key] = pd.concat([self._obs_values[key], values]).drop_duplicates()
        except KeyError:
            if isinstance(values, np.ndarray):
                self._obs_values[key] = pd.DataFrame(values, columns=[file_name]).T
            elif isinstance(values, pd.DataFrame):
                self._obs_values[key] = values


    def log_filename(self, frags):
        """
        """

        # Determine name
        path = os.path.normpath(frags.sam_file)
        file_name = path.split(os.sep)[-1]

        # Find all file names
        if self._anno is None:
            total_filenames = np.array([file_name])
        else:
            total_filenames = np.append(self._anno.index.values, file_name)
            self._anno = self._anno.reindex(total_filenames)

        # Record
        self.filenames.add(file_name)
        self.add_obs_values("length_dist", file_name, frags.length_dist())
        self.add_anno("n_fragments", file_name, len(frags.frags))

        if self.chrom_lengths is None:
            self.chrom_lengths = frags.genome
        if self.chroms is None:
            self.chroms = frags.chroms

        
