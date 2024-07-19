"""
Create matplotlib plots from cfDNA object
"""

import numpy as np
import matplotlib.pyplot as plt
import ngsfragments as ngs
import pandas as pd
from ..core.core import cfDNA

# Tell matplot lib to use 'editable text'
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def cnv_summary(cfdna_object: cfDNA,
                sample: str,
                add_wps: bool = False,
                show: bool = True,
                save: str = None) -> plt.figure:
    """
    Plot CNV summary

    Parameters
    ----------
        cfdna_object : cfDNA
            cfDNA object
        sample : str
            Sample name
        add_wps : bool
            Add WPS (default: False)
        show : bool
            Show plot (default: True)
        save : str
            Save plot to file (default: None)
    
    Returns
    -------
        fig : plt.figure
            matplotlib.pyplot.figure
    """

    # Create figure
    fig = ngs.plot.cnv_summary(cfdna_object,
                            sample,
                            add_wps = add_wps,
                            show = show,
                            save = save)

    return fig