"""InVEST Forage module.

InVEST Forage model developed from this design doc:
https://docs.google.com/document/d/10oJo43buEdJkFTZ0wYaW00EagSzs1oM7g_lBUc8URMI/edit#
"""
import os
import logging

LOGGER = logging.getLogger('natcap.invest.forage')

# we only have these types of soils
SOIL_TYPE_LIST = ['clay', 'silt', 'sand']

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
        args['clay_percent_path'] (string): path to raster representing per-pixel
            percent of soil component that is clay
        args['silt_percent_path'] (string): path to raster representing per-pixel
            percent of soil component that is silt
        args['sand_percent_path'] (string): path to raster representing per-pixel
            percent of soil component that is sand
        args['monthly_precip_path_pattern'] (string): path to the monthly
            precipitation path pattern. where the string <month> and <year>
            can be replaced with the number 1..12 for the month and integer
            year.  Example: if this value is given as:
            `./precip_dir/chirps-v2.0.<year>.<month>.tif, `starting_year` as
            2016, `starting_month` as 5, and `n_months` is 29, the model will
            expect to find files named
                 "./precip_dir/chirps-v2.0.2016.05.tif" to
                 "./precip_dir/chirps-v2.0.2018.09.tif"

    Returns:
        None.
    """
    LOGGER.info("model execute: %s", args)
    starting_month = int(args['starting_month'])
    starting_year = int(args['starting_year'])
    # this will be used to look up precip path given the model's month
    # timestep index that starts at 0
    precip_path_list = []
    # we'll use this to report any missing paths
    missing_path_list = []
    # build the list of
    for month_index in xrange(int(args['n_months'])):
        month_i = (starting_month + month_index - 1) % 12 + 1
        year = starting_year +  (starting_month + month_index - 1) // 12
        precip_path = args['monthly_precip_path_pattern'].replace(
            '<year>', str(year)).replace('<month>', '%.2d' % month_i)
        precip_path_list.append(precip_path)
        if not os.path.exists(precip_path):
            missing_path_list.append(precip_path)
    if missing_path_list:
        raise ValueError(
            "Couldn't find the following precipitation paths given the " +
            "pattern: %s\n\t" % args['monthly_precip_path_pattern'] +
            "\n\t".join(missing_path_list))

    # lookup to provide path to soil percent given soil type
    soil_path_type_map = {}
    for soil_type in SOIL_TYPE_LIST:
        soil_path_type_map[soil_type] = args['%s_percent_path' % soil_type]
        if not os.path.exists(soil_path_type_map[soil_type]):
            raise ValueError(
                "Couldn't find %s for %s" % (
                    soil_path_type_map[soil_type], soil_type))

