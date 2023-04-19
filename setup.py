from __future__ import division, print_function, absolute_import

# Python 3
MyFileNotFoundError = FileNotFoundError

# Library name
libname = "cfdna"

# Descriptions of package
SHORTDESC = "Python package for fragment manipulation for cfDNA"
#DESC = """A python package for fragment manipulation for cfDNA."""
with open("README.md", "r") as fh:
    long_description = fh.read()

# Directories (relative to the top-level directory where setup.py resides) in which to look for data files.
datadirs  = ("tests", "data")
# File extensions to be considered as data files. (Literal, no wildcards.)
dataexts  = (".py", ".pyx", ".pxd", ".c", ".h", ".h5", ".txt", ".bed")
# Standard documentation to detect (and package if it exists).
standard_docs     = ["README", "LICENSE", "TODO", "CHANGELOG", "AUTHORS"]
standard_doc_exts = [".md", ".rst", ".txt", ""]

#########################################################
# Init
#########################################################

# check for Python 2.7 or later
import sys
if sys.version_info < (2,7):
    sys.exit('Sorry, Python < 2.7 is not supported')

import os

from setuptools import setup, find_namespace_packages, find_packages


# Gather user-defined data files
datafiles = []
getext = lambda filename: os.path.splitext(filename)[1]
for datadir in datadirs:
    datafiles.extend( [(root, [os.path.join(root, f) for f in files if getext(f) in dataexts])
                       for root, dirs, files in os.walk(datadir)] )

# Add standard documentation (README et al.), if any, to data files
detected_docs = []
for docname in standard_docs:
    for ext in standard_doc_exts:
        filename = "".join( (docname, ext) )  # relative to the directory in which setup.py resides
        if os.path.isfile(filename):
            detected_docs.append(filename)
datafiles.append( ('.', detected_docs) )


# Extract __version__ from the package __init__.py
import ast
init_py_path = os.path.join(libname, '__init__.py')
version = '0.0.0'
try:
    with open(init_py_path) as f:
        for line in f:
            if line.startswith('__version__'):
                version = ast.parse(line).body[0].value.s
                break
        else:
            print( "WARNING: Version information not found in '%s', using placeholder '%s'" % (init_py_path, version), file=sys.stderr )
except MyFileNotFoundError:
    print( "WARNING: Could not find file '%s', using placeholder version information '%s'" % (init_py_path, version), file=sys.stderr )


#########################################################
# Call setup()
#########################################################

setup(
    name = "cfdna",
    version = version,
    author = "Kyle S. Smith",
    author_email = "kyle.smith@stjude.org",
    url = "https://github.com/kylessmtih/cfdna",
    description = SHORTDESC,
    long_description = long_description,
    long_description_content_type = "text/markdown",
    # CHANGE THIS
    license = "GPL2",
    # free-form text field
    platforms = ["Linux"],
    classifiers = [ "Development Status :: 4 - Beta",
                    "Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "Operating System :: POSIX :: Linux",
                    "Programming Language :: Python",
                    "Programming Language :: Python :: 3",
                    "Programming Language :: Python :: 3.4",
                    "Topic :: Scientific/Engineering",
                    "Topic :: Scientific/Engineering :: Mathematics",
                    "Topic :: Software Development :: Libraries",
                    "Topic :: Software Development :: Libraries :: Python Modules",
                    "Topic :: Scientific/Engineering :: Bio-Informatics"
                  ],
    setup_requires = ["numpy","ngsfragments","intervalframe","matplotlib","seaborn","bokeh","hmmCNV"],
    install_requires = ["numpy","ngsfragments","h5py","pandas"],
    provides = ["cfdna"],
    keywords = ["next generation sequencing fragment cfDNA cell free DNA"],
    packages = find_packages(),
    package_data={'cfdna': ['*.pxd', '*.pyx', '*.c', '*.h']},
    # Disable zip_safe
    zip_safe = False,
    # Custom data files not inside a Python package
    data_files = datafiles
)