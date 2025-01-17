from __future__ import absolute_import
from .core import cfDNA
from . import plot
from . import io
from . import commandline as cmd
from . import tools

from ngsfragments import Fragments


# This is extracted automatically by the top-level setup.py.
__version__ = '2.3.2'

__author__ = "Kyle S. Smith"

__doc__ = """\

API
======

Central class
-------------

.. autosummary::
   :toctree: .
   
   cfDNA


Plotting functions
------------------
    
.. autosummary::
   :toctree: .
   
   plot.plot_plt


IO functions
------------
    
.. autosummary::
   :toctree: .
   
   io.write_h5
   io.write_anno
   io.write_seg
   io.write_anno_seg
   io.read_sam
   io.read_h5


Tools
-----
    
.. autosummary::
   :toctree: .
   
   tools.nucleosome
   tools.cnv

"""