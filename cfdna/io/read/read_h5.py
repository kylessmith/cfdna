import pandas as pd
import numpy as np
import h5py
from projectframe import ProjectFrame

# Local imports
from ...core import cfDNA


def read_h5(h5_filename: str) -> cfDNA:
    """
    Read project frame from HDF5 file

    Parameters
    ----------
        h5_filename : str
            Name of HDF5 file to read from

    Returns
    -------
        cfdna_object : :class: `cfdna.cfDNA`
            cfDNA object
    """

    # Read project frame from HDF5 file
    cfdna_object = ProjectFrame.read_h5(h5_filename)

    return cfdna_object