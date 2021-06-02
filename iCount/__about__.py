"""Central place for package metadata."""

# NOTE: We use __title__ instead of simply __name__ since the latter would
#       interfere with a global variable __name__ denoting object's name.
__title__ = 'iCount-Mini'
__summary__ = 'Computational pipeline for analysis of iCLIP data'
__url__ = 'https://github.com/ulelab/iCount-Mini'

# Semantic versioning is used. For more information see:
# https://packaging.python.org/en/latest/distributing/#semantic-versioning-preferred
__version__ = '2.0.3'

__author__ = 'Ule Group, The Francis Crick Institute, London'
__email__ = 'jernej.ule@crick.ac.uk'

__license__ = 'MIT'
__copyright__ = '2016-2017, University of Ljubljana, Bioinformatics Laboratory'

__all__ = (
    '__title__', '__summary__', '__url__', '__version__', '__author__',
    '__email__', '__license__', '__copyright__',
)
