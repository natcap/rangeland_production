"""init module for rangeland_production."""

import os
import sys
import logging

import pkg_resources
import pygeoprocessing


LOGGER = logging.getLogger('rangeland_production')
LOGGER.addHandler(logging.NullHandler())
__all__ = ['local_dir',]

try:
    __version__ = pkg_resources.get_distribution(__name__).version
except pkg_resources.DistributionNotFound:
    # package is not installed.  Log the exception for debugging.
    LOGGER.exception('Could not load rangeland_production version information')
