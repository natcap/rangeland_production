"""InVEST Forage module.

InVEST Forage model developed from this design doc:
https://docs.google.com/document/d/10oJo43buEdJkFTZ0wYaW00EagSzs1oM7g_lBUc8URMI/edit#
"""
from __future__ import absolute_import

import os
import logging

from osgeo import gdal
from osgeo import ogr
import numpy

import pygeoprocessing
from . import utils
from . import validation

LOGGER = logging.getLogger('natcap.invest.forage')


def execute(args):
    """InVEST Forage Model.

    [model description]

    Parameters:
        args['workspace_dir'] (string): path to target output workspace.
        args['starting_month'] (int): what month to start reporting where
            the range 1..12 is equivalent to Jan..Dec.
        args['starting_year'] (int): what year to start runs. this value is
            used to notate outputs in the form [month_int]_[year]
        args['n_months'] (int): number of months to run model, the model run
            will start reporting in `args['starting_month']`.
        args['aoi_path'] (string): path to polygon vector indicating the
            desired spatial extent of the model. This has the effect of
            clipping the computational area of the input datasets to be the
            area intersected by this polygon.

    Returns:
        None.
    """
    LOGGER.info("model execute: %s", args)
