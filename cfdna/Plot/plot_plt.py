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
                key: str,
                show: bool = True,
                save: bool = None) -> plt.figure:
    """
    Plot summary metrics

    Parameters
    ----------
        cfdna_object : cfDNA
            cfDNA object
        key : str
            Sample name
        show : bool
            Show plot
        save : str
            Save plot
    
    Returns
    -------
        fig
            matplotlib.pyplot.figure
    """

    # Create figure
    fig = ngs.plot.plot_plt.cnv_summary(cfdna_object,
                                        key,
                                        show = show,
                                        save = save)

    return fig