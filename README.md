# Cell-free DNA analysis toolkit

[![Build Status](https://travis-ci.org/kylessmith/cfdna.svg?branch=master)](https://travis-ci.org/kylessmith/cfdna) [![PyPI version](https://badge.fury.io/py/cfdna.svg)](https://badge.fury.io/py/cfdna)

This is a Python package for easy and efficient cell-free
DNA analysis.

All citations should reference to [original paper][paper].


## Install

If you dont already have numpy and scipy installed, it is best to download
`Anaconda`, a python distribution that has them included.  
```
    https://continuum.io/downloads
```

Create environment and pip install:
```
    conda create -n cfdna python=3.11.3
	conda activate cfdna
	pip install cfdna==2.0.5
```

Dependencies can be installed by:

```
    pip install -r requirements.txt
```

PyPI install, presuming you have all its requirements installed:
```
    pip install cfdna
```

## Usage

```python
import cfdna as cf
import numpy as np

# Create data
cfdata = cf.cfDNA()
frags = cf.io.readBAM("test.bam")

# Call CNVs
cf.proc.call_cnv_pipline(cfdata, frags)

# Plot CNVs
cf.pl.plot_cnv(cfdata, "test.bam")
```

Run from the commadline:

```
    python -m cfdna callCNVs --bam test.bam --segs --genome hg19 --bin_size 100000
```

This will output a .png plot and seg file.

[paper]: https://www.cell.com/cancer-cell/pdfExtended/S1535-6108(21)00501-8