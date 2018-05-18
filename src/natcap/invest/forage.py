"""InVEST Forage module.

InVEST Forage model developed from this design doc:
https://docs.google.com/document/d/10oJo43buEdJkFTZ0wYaW00EagSzs1oM7g_lBUc8URMI/edit#
"""
import os
import logging
import tempfile
import shutil

from osgeo import ogr
from osgeo import gdal
import re
import numpy

import pygeoprocessing
from natcap.invest import utils

LOGGER = logging.getLogger('natcap.invest.forage')

# we only have these types of soils
SOIL_TYPE_LIST = ['clay', 'silt', 'sand']

# state variables and parameters take their names from Century
# _SITE_STATE_VARIABLE_FILES contains state variables that are a
# property of the site, including:
# carbon in each soil compartment
# (structural, metabolic, som1, som2, som3) and layer (1=surface, 2=soil)
# e.g., som2c_2 = carbon in soil som2;
# N and P in each soil layer and compartment (1=N, 2=P)
# e.g., som2e_1_1 = N in surface som2, som2e_1_2 = P in surface som2;
# water in each soil layer, asmos_<layer>
# state variables fully described in this table:
# https://docs.google.com/spreadsheets/d/1TGCDOJS4nNsJpzTWdiWed390NmbhQFB2uUoMs9oTTYo/edit?usp=sharing
_SITE_STATE_VARIABLE_FILES = {
    'metabc_1_path': 'metabc_1.tif',
    'metabc_2_path': 'metabc_2.tif',
    'som1c_1_path': 'som1c_1.tif',
    'som1c_2_path': 'som1c_2.tif',
    'som2c_1_path': 'som2c_1.tif',
    'som2c_2_path': 'som2c_2.tif',
    'som3c_path': 'som3c.tif',
    'strucc_1_path': 'strucc_1.tif',
    'strucc_2_path': 'strucc_2.tif',
    'metabe_1_1_path': 'metabe_1_1.tif',
    'metabe_2_1_path': 'metabe_2_1.tif',
    'minerl_1_1_path': 'minerl_1_1.tif',
    'som1e_1_1_path': 'som1e_1_1.tif',
    'som1e_2_1_path': 'som1e_2_1.tif',
    'som2e_1_1_path': 'som2e_1_1.tif',
    'som2e_2_1_path': 'som2e_2_1.tif',
    'som3e_1_path': 'som3e_1.tif',
    'struce_1_1_path': 'struce_1_1.tif',
    'struce_2_1_path': 'struce_2_1.tif',
    'metabe_1_2_path': 'metabe_1_2.tif',
    'metabe_2_2_path': 'metabe_2_2.tif',
    'plabil_path': 'plabil.tif',
    'secndy_2_path': 'secndy_2.tif',
    'som1e_1_2_path': 'som1e_1_2.tif',
    'som1e_2_2_path': 'som1e_2_2.tif',
    'som2e_1_2_path': 'som2e_1_2.tif',
    'som2e_2_2_path': 'som2e_2_2.tif',
    'som3e_1_path': 'som3e_1.tif',
    'struce_1_2_path': 'struce_1_2.tif',
    'struce_2_2_path': 'struce_2_2.tif',
    'asmos_1_path': 'asmos_1.tif',
    'asmos_2_path': 'asmos_2.tif',
    'asmos_3_path': 'asmos_3.tif',
    'asmos_4_path': 'asmos_4.tif',
    'asmos_5_path': 'asmos_5.tif',
    'asmos_6_path': 'asmos_6.tif',
    'asmos_7_path': 'asmos_7.tif',
    'asmos_8_path': 'asmos_8.tif',
    'asmos_9_path': 'asmos_9.tif',
    }

# _PFT_STATE_VARIABLES contains state variables that are a
# property of a PFT, including:
# carbon, nitrogen, and phosphorous in aboveground biomass
# where 1=N, 2=P
# e.g. aglivc = C in aboveground live biomass,
# aglive_1 = N in aboveground live biomass;
# carbon, nitrogen, and phosphorous in aboveground standing dead
# biomass, stdedc and stdede;
# carbon, nitrogen and phosphorous in belowground live biomass,
# aglivc and aglive
# state variables fully described in this table:
# https://docs.google.com/spreadsheets/d/1TGCDOJS4nNsJpzTWdiWed390NmbhQFB2uUoMs9oTTYo/edit?usp=sharing
_PFT_STATE_VARIABLES = [
    'aglivc', 'bglivc', 'stdedc', 'aglive_1', 'bglive_1',
    'stdede_1', 'aglive_2', 'bglive_2', 'stdede_2',
    ]

# intermediate parameters that do not change between timesteps,
# including field capacity and wilting point of each soil layer,
# coefficients describing effect of soil texture on decomposition
# rates

_PERSISTENT_PARAMS_FILES = {
    'afiel_1_path': 'afiel_1.tif',
    'afiel_2_path': 'afiel_2.tif',
    'afiel_3_path': 'afiel_3.tif',
    'afiel_4_path': 'afiel_4.tif',
    'afiel_5_path': 'afiel_5.tif',
    'afiel_6_path': 'afiel_6.tif',
    'afiel_7_path': 'afiel_7.tif',
    'afiel_8_path': 'afiel_8.tif',
    'afiel_9_path': 'afiel_9.tif',
    'awilt_1_path': 'awilt_1.tif',
    'awilt_2_path': 'awilt_2.tif',
    'awilt_3_path': 'awilt_3.tif',
    'awilt_4_path': 'awilt_4.tif',
    'awilt_5_path': 'awilt_5.tif',
    'awilt_6_path': 'awilt_6.tif',
    'awilt_7_path': 'awilt_7.tif',
    'awilt_8_path': 'awilt_8.tif',
    'awilt_9_path': 'awilt_9.tif',
    'wc_path': 'wc.tif',
    'eftext_path': 'eftext.tif',
    'p1co2_2_path': 'p1co2_2.tif',
    'fps1s3_path': 'fps1s3.tif',
    'orglch_path': 'orglch.tif',
    'fps2s3_path': 'fps2s3.tif',
    'rnewas_1_1_path': 'rnewas_1_1.tif',
    'rnewas_2_1_path': 'rnewas_2_1.tif',
    'rnewas_1_2_path': 'rnewas_1_2.tif',
    'rnewas_2_2_path': 'rnewas_2_2.tif',
    'rnewbs_1_1_path': 'rnewbs_1_1.tif',
    'rnewbs_1_2_path': 'rnewbs_1_2.tif',
    'rnewbs_2_1_path': 'rnewbs_2_1.tif',
    'rnewbs_2_2_path': 'rnewbs_2_2.tif',
    }
# values that are updated once per year
_YEARLY_FILES = {
    'annual_precip_path': 'annual_precip.tif',
    'baseNdep_path': 'baseNdep.tif',
    }

# Target nodata is for general rasters that are positive, and _IC_NODATA are
# for rasters that are any range
_TARGET_NODATA = -1.0
_IC_NODATA = numpy.finfo('float32').min


def execute(args):
    """InVEST Forage Model.

    [model description]

    Parameters:
        args['workspace_dir'] (string): path to target output workspace.
        args['results_suffix'] (string): (optional) string to append to any
            output file names
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
        args['bulk_density_path'] (string): path to bulk density raster.
        args['ph_path'] (string): path to soil pH raster.
        args['clay_proportion_path'] (string): path to raster representing
            per-pixel proportion of soil component that is clay
        args['silt_proportion_path'] (string): path to raster representing
            per-pixel proportion of soil component that is silt
        args['sand_proportion_path'] (string): path to raster representing
            per-pixel proportion of soil component that is sand
        args['monthly_precip_path_pattern'] (string): path to the monthly
            precipitation path pattern. where the string <month> and <year>
            can be replaced with the number 1..12 for the month and integer
            year. The model requires at least 12 months of precipitation and
            expects to find a precipitation file input for
            every month of the simulation, so the number of precipitation
            files should be the maximum of 12 and `n_months`. The <month>
            value in input files must be two digits.
            Example: if this value is given as:
            `./precip_dir/chirps-v2.0.<year>.<month>.tif, `starting_year` as
            2016, `starting_month` as 5, and `n_months` is 29, the model will
            expect to find files named
                 "./precip_dir/chirps-v2.0.2016.05.tif" to
                 "./precip_dir/chirps-v2.0.2018.09.tif"
        args['min_temp_path_pattern'] (string): path to monthly
            temperature data pattern where <month> can be replaced with the
            number 1..12 when the simulation needs a monthly temperature
            input. The model expects to find only one minimum temperature
            file input for each month of the year, so the number of minimum
            temperature input files could be less than `n_months`. The
            <month> value in input files must be two digits.
            Example: if this value is given as
            `./temperature/min_temperature_<month>.tif', `starting_month` as
            5, and `n_months` is 29, the model will expect to find files
            named
                "./temperature/min_temperature_01.tif to
                "./temperature/min_temperature_12.tif"
        args['max_temp_path_pattern'] (string): path to monthly
            temperature data pattern where <month> can be replaced with the
            number 1..12 when the simulation needs a monthly temperature
            input. The model expects to find only one maximum temperature
            file input for each month of the year, so the number of maximum
            temperature input files could be less than `n_months`. The
            <month> value in input files must be two digits.
            Example: if this value is given as
            `./temperature/max_temperature_<month>.tif', `starting_month` as
            5, and `n_months` is 29, the model will expect to find files
            named
                "./temperature/max_temperature_01.tif to
                "./temperature/max_temperature_12.tif"
        args['site_param_table'] (string): path to csv file giving site
            parameters. This file must contain a column named "site" that
            contains unique integers. These integer values correspond to site
            type identifiers which are values in the site parameter spatial
            index raster. Other required fields for this table are site and
            "fixed" parameters from the Century model, i.e., the parameters
            in the Century input files site.100 and fix.100.
        args['site_param_spatial_index_path'] (string): path to a raster file
            that indexes site parameters, indicating which set of site
            parameter values should apply at each pixel in the raster. The
            raster should be composed of integers that correspond to values in
            the field "site" in `site_param_table`.
        args['veg_trait_path'] (string): path to csv file giving vegetation
            traits for each plant functional type available for grazing. This
            file must contain a column named "PFT" that contains unique
            integers. These integer values correspond to PFT identifiers of
            veg spatial composition rasters. Other required fields for this
            table are vegetation input parameters from the Century model, for
            example maximum intrinsic growth rate, optimum temperature for
            production, minimum C/N ratio, etc.
        args['veg_spatial_composition_path_pattern'] (string): path to
            vegetation rasters, one per plant functional type available for
            grazing, where <PFT> can be replaced with an integer that is
            indexed in the veg trait csv.
            Example: if this value is given as `./vegetation/pft_<PFT>.tif`
            and the directory `./vegetation/` contains these files:
                "pft_1.tif"
                "pft_12.tif"
                "pft_50.tif",
            then the "PFT" field in the vegetation trait table must contain
            the values 1, 12, and 50.
        args['animal_trait_path'] (string): path to csv file giving animal
            traits for each animal type - number - duration combination. This
            table must contain a column named "animal_id" that contains unique
            integers. These integer values correspond to features in the
            animal management layer.
            Other required fields in this table are:
                type (allowable values: B_indicus, B_taurus,
                    indicus_x_taurus, sheep, camelid, hindgut_fermenter)
                sex (allowable values: entire_m, castrate, lac_female,
                    nonlac_female)
                age (days)
                weight (kg)
                SRW (standard reference weight, kg; the weight of a mature
                    female in median condition)
                SFW (standard fleece weight, kg; the average weight of fleece
                    of a mature adult; for sheep only)
                birth_weight (kg)
                num_animals (number of animals grazing inside the polygon)
                grz_months (a string of integers, separated by ','; months of
                    the simulation when animals are present,
                    relative to `starting_month`. For example, if `n_months`
                    is 3, and animals are present during the entire simulation
                    period, `grz_months` should be "1,2,3")
        args['animal_mgmt_layer_path'] (string): path to animal vector inputs
            giving the location of grazing animals. Must have a field named
            "animal_id", containing unique integers that correspond to the
            values in the "animal_id" column of the animal trait csv.
        args['initial_conditions_dir'] (string): optional path to directory
            containing initial conditions. If supplied, this directory must
            contain a series of rasters with initial values for each PFT and
            for the site.
                Required rasters for each PFT:
                    initial variables that are a property of PFT in the table
                    https://docs.google.com/spreadsheets/d/1TGCDOJS4nNsJpzTWdiWed390NmbhQFB2uUoMs9oTTYo/edit?usp=sharing
                    e.g., aglivc_<PFT>.tif
                Required for the site:
                    initial variables that are a property of site in the table
                    https://docs.google.com/spreadsheets/d/1TGCDOJS4nNsJpzTWdiWed390NmbhQFB2uUoMs9oTTYo/edit?usp=sharing

    Returns:
        None.
    """
    LOGGER.info("model execute: %s", args)
    starting_month = int(args['starting_month'])
    starting_year = int(args['starting_year'])
    n_months = int(args['n_months'])
    # this will be used to look up precip path given the model's month
    # timestep index that starts at 0
    precip_path_list = []
    # this set will build up the integer months that are used so we can index
    # the mwith temperature later
    temperature_month_set = set()

    # this dict will be used to build the set of input rasters associated with
    # a reasonable lookup ID so we can have a nice dataset to align for raster
    # stack operations
    base_align_raster_path_id_map = {}
    # we'll use this to report any missing paths
    missing_precip_path_list = []
    for month_index in xrange(n_months):
        month_i = (starting_month + month_index - 1) % 12 + 1
        temperature_month_set.add(month_i)
        year = starting_year + (starting_month + month_index - 1) // 12
        precip_path = args[
            'monthly_precip_path_pattern'].replace(
                '<year>', str(year)).replace('<month>', '%.2d' % month_i)
        base_align_raster_path_id_map['precip_%d' % month_index] = precip_path
        precip_path_list.append(precip_path)
        if not os.path.exists(precip_path):
            missing_precip_path_list.append(precip_path)
    if missing_precip_path_list:
        raise ValueError(
            "Couldn't find the following precipitation paths given the "
            + "pattern: %s\n\t" % args['monthly_precip_path_pattern']
            + "\n\t".join(missing_precip_path_list))
    # the model requires 12 months of precipitation data to calculate
    # atmospheric N deposition and potential production from annual precip
    n_precip_months = int(args['n_months'])
    if n_precip_months < 12:
        month_index = int(args['n_months'])
        while month_index <= 12:
            month_i = (starting_month + month_index - 1) % 12 + 1
            year = starting_year + (starting_month + month_index - 1) // 12
            precip_path = args['monthly_precip_path_pattern'].replace(
                '<year>', str(year)).replace('<month>', '%.2d' % month_i)
            base_align_raster_path_id_map['precip_%d' % month_index] = precip_path
            precip_path_list.append(precip_path)
            if os.path.exists(precip_path):
                n_precip_months = n_precip_months + 1
            month_index = month_index + 1
    if n_precip_months < 12:
        raise ValueError("At least 12 months of precipitation data required")

    # this list will be used to record any expected files that are not found
    for substring in ['min', 'max']:
        missing_temperature_path_list = []
        for month_i in temperature_month_set:
            monthly_temp_path = args[
                '%s_temp_path_pattern' % substring].replace(
                    '<month>', '%.2d' % month_i)
            base_align_raster_path_id_map[
                '%s_temp_%d' % (substring, month_i)] = monthly_temp_path
            if not os.path.exists(monthly_temp_path):
                missing_min_temperature_path_list.append(monthly_temp_path)
        if missing_temperature_path_list:
            raise ValueError(
                "Couldn't find the following temperature raster paths"
                + " given the pattern: %s\n\t" % args[
                    '%s_temp_path_pattern' % substring]
                + "\n\t".join(missing_temperature_path_list))

    # lookup to provide path to soil percent given soil type
    for soil_type in SOIL_TYPE_LIST:
        base_align_raster_path_id_map[soil_type] = (
            args['%s_proportion_path' % soil_type])
        if not os.path.exists(base_align_raster_path_id_map[soil_type]):
            raise ValueError(
                "Couldn't find %s for %s" % (
                    base_align_raster_path_id_map[soil_type], soil_type))
    base_align_raster_path_id_map['bulk_d_path'] = args['bulk_density_path']
    base_align_raster_path_id_map['ph_path'] = args['ph_path']

    # make sure site parameters exist for each site type identifier
    base_align_raster_path_id_map['site_index'] = (
        args['site_param_spatial_index_path'])
    n_bands = pygeoprocessing.get_raster_info(
        args['site_param_spatial_index_path'])['n_bands']
    if n_bands > 1:
        raise ValueError(
            'Site spatial index raster must contain only one band')
    site_datatype = pygeoprocessing.get_raster_info(
        args['site_param_spatial_index_path'])['datatype']
    if site_datatype not in [1, 2, 3, 4, 5]:
        raise ValueError('Site spatial index raster must be integer type')

    # get unique values in site param raster
    site_index_set = set()
    for offset_map, raster_block in pygeoprocessing.iterblocks(
            args['site_param_spatial_index_path']):
        site_index_set.update(numpy.unique(raster_block))
    site_nodata = pygeoprocessing.get_raster_info(
        args['site_param_spatial_index_path'])['nodata'][0]
    if site_nodata in site_index_set:
        site_index_set.remove(site_nodata)
    site_param_table = utils.build_lookup_from_csv(
        args['site_param_table'], 'site')
    missing_site_index_list = list(
        site_index_set.difference(site_param_table.iterkeys()))
    if missing_site_index_list:
        raise ValueError(
            "Couldn't find parameter values for the following site "
            + "indices: %s\n\t" + ", ".join(missing_site_index_list))

    # make sure veg traits exist for each pft raster
    pft_dir = os.path.dirname(args['veg_spatial_composition_path_pattern'])
    pft_basename = os.path.basename(
        args['veg_spatial_composition_path_pattern'])
    files = [
        f for f in os.listdir(pft_dir) if os.path.isfile(
            os.path.join(pft_dir, f))]
    pft_regex = re.compile(pft_basename.replace('<PFT>', r'(\d+)'))
    pft_matches = [
        m for m in [pft_regex.search(f) for f in files] if m is not None]
    pft_id_set = set([int(m.group(1)) for m in pft_matches])
    for pft_i in pft_id_set:
        pft_path = args['veg_spatial_composition_path_pattern'].replace(
            '<PFT>', '%d' % pft_i)
        base_align_raster_path_id_map['pft_%d' % pft_i] = pft_path
    veg_trait_table = utils.build_lookup_from_csv(
        args['veg_trait_path'], 'PFT')
    missing_pft_trait_list = pft_id_set.difference(veg_trait_table.iterkeys())
    if missing_pft_trait_list:
        raise ValueError(
            "Couldn't find trait values for the following plant functional "
            + "types: %s\n\t" + ", ".join(missing_pft_trait_list))

    # track separate state variable files for each PFT
    pft_sv_dict = {}
    for pft_index in pft_id_set:
        for sv in _PFT_STATE_VARIABLES:
            pft_sv_dict['{}_{}_path'.format(
                sv, pft_index)] = '{}_{}.tif'.format(sv, pft_index)

    # make sure animal traits exist for each feature in animal management
    # layer
    anim_id_list = []
    driver = ogr.GetDriverByName('ESRI Shapefile')
    datasource = driver.Open(args['animal_mgmt_layer_path'], 0)
    layer = datasource.GetLayer()
    for feature in layer:
        anim_id_list.append(feature.GetField('animal_id'))

    anim_trait_table = utils.build_lookup_from_csv(
        args['animal_trait_path'], 'animal_id')
    missing_animal_trait_list = set(
        anim_id_list).difference(anim_trait_table.iterkeys())
    if missing_animal_trait_list:
        raise ValueError(
            "Couldn't find trait values for the following animal "
            + "ids: %s\n\t" + ", ".join(missing_animal_trait_list))

    # find the smallest target_pixel_size
    target_pixel_size = min(*[
        pygeoprocessing.get_raster_info(path)['pixel_size']
        for path in base_align_raster_path_id_map.values()],
        key=lambda x: (abs(x[0]), abs(x[1])))
    LOGGER.info(
        "smallest pixel size of all raster inputs: %s", target_pixel_size)

    # set up a dictionary that uses the same keys as
    # 'base_align_raster_path_id_map' to point to the clipped/resampled
    # rasters to be used in raster calculations for the model.
    aligned_raster_dir = os.path.join(
        args['workspace_dir'], 'aligned_inputs')
    aligned_inputs = dict([(key, os.path.join(
        aligned_raster_dir, 'aligned_%s' % os.path.basename(path)))
        for key, path in base_align_raster_path_id_map.iteritems()])

    # align all the base inputs to be the minimum known pixel size and to
    # only extend over their combined intersections
    source_input_path_list = [
        base_align_raster_path_id_map[k] for k in sorted(
            base_align_raster_path_id_map.iterkeys())]
    aligned_input_path_list = [
        aligned_inputs[k] for k in sorted(aligned_inputs.iterkeys())]

    file_suffix = utils.make_suffix_string(args, 'results_suffix')

    # if initial conditions are supplied, use them to initialize all state
    # variables
    if args['initial_conditions_dir']:
        LOGGER.info(
            "setting initial conditions from this directory: %s",
            args['initial_conditions_dir'])

        # check that all necessary state variables are supplied
        missing_initial_values = []
        # align initial state variables to resampled inputs
        resample_initial_path_map = {}
        for sv in _SITE_STATE_VARIABLE_FILES.iterkeys():
            sv_path = os.path.join(
                args['initial_conditions_dir'],
                _SITE_STATE_VARIABLE_FILES[sv])
            resample_initial_path_map[sv] = sv_path
            if not os.path.exists(sv_path):
                missing_initial_values.append(sv_path)
        for pft_index in pft_id_set:
            for sv in _PFT_STATE_VARIABLES:
                sv_key = '{}_{}_path'.format(sv, pft_index)
                sv_path = os.path.join(
                    args['initial_conditions_dir'],
                    '{}_{}.tif'.format(sv, pft_index))
                resample_initial_path_map[sv_key] = sv_path
                if not os.path.exists(sv_path):
                    missing_initial_values.append(sv_path)
        if missing_initial_values:
            raise ValueError(
                "Couldn't find the following required initial values: "
                + "\n\t".join(missing_initial_values))

        # align and resample initialization rasters with inputs
        sv_dir = os.path.join(args['workspace_dir'], 'state_variables_m-1')
        aligned_initial_path_map = dict(
            [(key, os.path.join(sv_dir, os.path.basename(path)))
                for key, path in resample_initial_path_map.iteritems()])
        initial_path_list = (
            [resample_initial_path_map[k] for k in
                sorted(resample_initial_path_map.iterkeys())]
            + source_input_path_list)
        aligned_initial_path_list = (
            [aligned_initial_path_map[k] for k in
                sorted(aligned_initial_path_map.iterkeys())]
            + aligned_input_path_list)
        pygeoprocessing.align_and_resize_raster_stack(
            initial_path_list, aligned_initial_path_list,
            ['nearest'] * len(initial_path_list),
            target_pixel_size, 'intersection', raster_align_index=0)
        sv_reg = aligned_initial_path_map
    else:
        LOGGER.info("initial conditions not supplied")
        LOGGER.info("aligning base raster inputs")
        pygeoprocessing.align_and_resize_raster_stack(
            source_input_path_list, aligned_input_path_list,
            ['nearest'] * len(aligned_inputs),
            target_pixel_size, 'intersection',
            base_vector_path_list=[args['aoi_path']])
        # TODO add spin-up or initialization with Burke's equations
        raise ValueError("Initial conditions must be supplied")

    # Initialization
    # calculate persistent intermediate parameters that do not change during
    # the simulation
    persist_param_dir = os.path.join(
        args['workspace_dir'], 'intermediate_parameters')
    utils.make_directories([persist_param_dir])
    pp_reg = utils.build_file_registry(
        [(_PERSISTENT_PARAMS_FILES, persist_param_dir)], file_suffix)

    # calculate field capacity and wilting point
    LOGGER.info("Calculating field capacity and wilting point")
    _afiel_awilt(
        aligned_inputs['site_index'], site_param_table,
        sv_reg['som1c_2_path'], sv_reg['som2c_2_path'], sv_reg['som3c_path'],
        aligned_inputs['sand'], aligned_inputs['silt'],
        aligned_inputs['clay'], aligned_inputs['bulk_d_path'], pp_reg)

    # calculate other persistent parameters
    LOGGER.info("Calculating persistent parameters")
    _persistent_params(
        aligned_inputs['site_index'], site_param_table,
        aligned_inputs['sand'], aligned_inputs['clay'], pp_reg)
<<<<<<< working copy

    # calculate required ratios for decomposition of structural material
    LOGGER.info("Calculating required ratios for structural decomposition")
    _structural_ratios(aligned_inputs['site_index'], site_param_table,
        sv_reg, pp_reg)

    # make yearly directory for values that are updated every twelve months
    year_dir = tempfile.mkdtemp()
    year_reg = dict([(key, os.path.join(year_dir, path)) for key, path in
        _YEARLY_FILES.iteritems()])

    # make monthly directory for monthly intermediate parameters that are
    # shared between submodels, but do not need to be saved as output
    month_temp_dir = tempfile.mkdtemp()
    month_reg = dict([(key, os.path.join(month_temp_dir, path)) for
        key, path in _SITE_INTERMEDIATE_PARAMS.iteritems()])
    for pft_index in pft_id_set:
        for val in _PFT_INTERMEDIATE_PARAMS:
            month_reg['{}_{}_path'.format(
                val, pft_index)] = '{}_{}.tif'.format(val, pft_index)

    # Main simulation loop
    # for each step in the simulation
    for month_index in xrange(n_months):
        if (month_index % 12) == 0:
            # Update yearly quantities
            _yearly_tasks(aligned_inputs['site_index'], site_param_table,
                aligned_inputs, month_index, year_reg)

        month_i = (starting_month + month_index - 1) % 12 + 1
        year = starting_year + (starting_month + month_index - 1) // 12

        # make new folders for state variables during this step
        sv_dir = os.path.join(
            args['workspace_dir'], 'state_variables_m%d' % month_index)
        utils.make_directories([sv_dir])

        # track state variables from previous step
        prev_sv_reg = sv_reg
        sv_reg = utils.build_file_registry(
            [(_SITE_STATE_VARIABLE_FILES, sv_dir),
                (pft_sv_dict, sv_dir)], file_suffix)

        # update state variables from previous month
        LOGGER.info(
            "Main simulation loop: month %d of %d" % (
                month_index, n_months))


def _afiel_awilt(site_index_path, site_param_table, som1c_2_path,
                 som2c_2_path, som3c_path, sand_path, silt_path, clay_path,
                 bulk_d_path, pp_reg):
    """Calculate field capacity and wilting point for each soil layer.

    Computations based on Gupta and Larson 1979, 'Estimating soil and water
    retention characteristics from particle size distribution, organic
    matter percent and bulk density'. Water Resources Research 15:1633.
    Field capacity is calculated for -0.33 bar; wilting point is
    calculated for water content at -15 bars.

    Parameters:
        site_index_path (string): path to site spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters including 'edepth' field
        som1c_2_path (string): path to the state variable 'som1c_2',
            active organic soil carbon
        som2c_2_path (string): path to the state variable 'som2c_2',
            slow organic soil carbon
        som3c_path (string): path to the state variable 'som3c',
            passive organic soil carbon
        sand_path (string): path to raster containing proportion sand in soil
        silt_path (string): path to raster containing proportion silt in soil
        clay_path (string): path to raster containing proportion clay in soil
        bulk_d_path (string): path to raster containing bulk density of soil
        pp_reg (dict): map of key, path pairs giving paths to persistent
            intermediate parameters that do not change over the course of
            the simulation.

    Modifies the rasters pp_reg['afiel_<layer>'] and pp_reg['awilt_<layer>']
        for all soil layers.

    Returns:
        None
    """
    # temporary intermediate rasters for this calculation
    temp_dir = tempfile.mkdtemp()
    edepth_path = os.path.join(temp_dir, 'edepth.tif')
    ompc_path = os.path.join(temp_dir, 'ompc.tif')

    def calc_ompc(som1c_2, som2c_2, som3c, bulkd, edepth):
        """Estimate total soil organic matter.

        Parameters:
            som1c_2 (numpy.ndarray): active organic soil carbon
            som2c_2 (numpy.ndarray): slow organic soil carbon
            som3c (numpy.ndarray): passive organic soil carbon
            bulkd (numpy.ndarray): bulk density of soil
            edepth (numpy.ndarray): parameter, depth of soil for this
                calculation
        From line 222, Prelim.f
        Returns:
            ompc, total soil organic matter weighted by bulk
            density.
        """
        ompc = numpy.empty(som1c_2.shape, dtype=numpy.float32)
        ompc[:] = _TARGET_NODATA
        valid_mask = (
            (som1c_2 > 0) & (som2c_2 > 0) & (som3c > 0) & (bulkd > 0))
        ompc[valid_mask] = (
            (som1c_2[valid_mask] + som2c_2[valid_mask]
                + som3c[valid_mask]) * 1.724
            / (10000. * bulkd[valid_mask] * edepth[valid_mask]))
        return ompc

    def calc_afiel(sand, silt, clay, ompc, bulkd):
        """Calculate field capacity for one soil layer.

        Field capacity, maximum soil moisture retention capacity,
        from Gupta and Larson 1979, 'Estimating soil and water
        retention characteristics from particle size distribution,
        organic matter percent and bulk density'. Water Resources
        Research 15:1633.

        Parameters:
            sand (numpy.ndarray): proportion sand in soil
            silt (numpy.ndarray): proportion silt in soil
            clay (numpy.ndarray): proportion clay in soil
            ompc (numpy.ndarray): estimated total soil organic matter
            bulkd (numpy.ndarray): bulk density of soil

        Returns:
            afiel, field capacity for this soil layer
        """
        afiel = numpy.empty(sand.shape, dtype=numpy.float32)
        afiel[:] = _TARGET_NODATA
        valid_mask = ompc != _TARGET_NODATA
        afiel[valid_mask] = (
            0.3075 * sand[valid_mask] + 0.5886 * silt[valid_mask]
            + 0.8039 * clay[valid_mask] + 2.208E-03 * ompc[valid_mask]
            + -0.1434 * bulkd[valid_mask])
        return afiel

    def calc_awilt(sand, silt, clay, ompc, bulkd):
        """Calculate wilting point for one soil layer.

        Wilting point, minimum soil water required by plants before
        wilting, from Gupta and Larson 1979, 'Estimating soil and
        water retention characteristics from particle size distribution,
        organic matter percent and bulk density'. Water Resources
        Research 15:1633.

        Parameters:
            sand (numpy.ndarray): proportion sand in soil
            silt (numpy.ndarray): proportion silt in soil
            clay (numpy.ndarray): proportion clay in soil
            ompc (numpy.ndarray): estimated total soil organic matter
            bulkd (numpy.ndarray): bulk density of soil

        Returns:
            awilt, wilting point for this soil layer
        """
        awilt = numpy.empty(sand.shape, dtype=numpy.float32)
        awilt[:] = _TARGET_NODATA
        valid_mask = ompc != _TARGET_NODATA
        awilt[valid_mask] = (
            -0.0059 * sand[valid_mask] + 0.1142 * silt[valid_mask]
            + 0.5766 * clay[valid_mask] + 2.228E-03 * ompc[valid_mask]
            + 0.02671 * bulkd[valid_mask])
        return awilt

    def decrement_ompc(ompc):
        """Decrease estimated organic matter in each subsequent layer."""
        return ompc * 0.85

    # temporary raster from site parameter 'edepth'
    site_to_edepth = dict(
        [(site_code, float(table['edepth'])) for
         (site_code, table) in site_param_table.iteritems()])

    pygeoprocessing.reclassify_raster(
        (site_index_path, 1), site_to_edepth, edepth_path, gdal.GDT_Float32,
        _TARGET_NODATA)

    # estimate total soil organic matter
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            som1c_2_path, som2c_2_path, som3c_path,
            bulk_d_path, edepth_path]],
        calc_ompc, ompc_path, gdal.GDT_Float32, _TARGET_NODATA)

    # calculate field capacity and wilting point for each soil layer,
    # decreasing organic matter content by 85% with each layer
    for lyr in xrange(1, 10):
        afiel_path = pp_reg['afiel_{}_path'.format(lyr)]
        awilt_path = pp_reg['awilt_{}_path'.format(lyr)]
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sand_path, silt_path, clay_path, ompc_path, bulk_d_path]],
            calc_afiel, afiel_path, gdal.GDT_Float32, _TARGET_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sand_path, silt_path, clay_path, ompc_path, bulk_d_path]],
            calc_awilt, awilt_path, gdal.GDT_Float32, _TARGET_NODATA)
        ompc_dec_path = os.path.join(temp_dir, 'ompc{}.tif'.format(lyr))
        pygeoprocessing.raster_calculator(
            [(ompc_path, 1)], decrement_ompc, ompc_dec_path, gdal.GDT_Float32,
            _TARGET_NODATA)
        ompc_path = ompc_dec_path

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _persistent_params(site_index_path, site_param_table, sand_path,
                       clay_path, pp_reg):
    """Calculate persistent parameters.

    The calculated values do not change over the course of the simulation.

    Parameters:
        site_index_path (string): path to site spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        sand_path (string): path to raster containing proportion sand in soil
        clay_path (string): path to raster containing proportion clay in soil
        pp_reg (dict): map of key, path pairs giving paths to persistent
            intermediate parameters that do not change over the course of
            the simulation.

    Modifies the persistent parameter rasters indexed by the following
    keys:
        pp_reg['eftext_path']
        pp_reg['p1co2_2_path']
        pp_reg['fps1s3_path']
        pp_reg['fps2s3_path']

    Returns:
        None
    """
    # temporary intermediate rasters for these calculations
    temp_dir = tempfile.mkdtemp()
    param_val_dict = {}
    for val in[
            'peftxa', 'peftxb', 'p1co2a_2', 'p1co2b_2', 'ps1s3_1',
            'ps1s3_2', 'ps2s3_1', 'ps2s3_2', 'omlech_1', 'omlech_2']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for (
                site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_val, target_path, gdal.GDT_Float32,
            _TARGET_NODATA)

    sand_nodata = pygeoprocessing.get_raster_info(
        sand_path)['nodata'][0]
    clay_nodata = pygeoprocessing.get_raster_info(
        clay_path)['nodata'][0]

    def calc_wc(afiel_1, awilt_1):
        """Calculate water content of soil layer 1, which is used to predict
        potential production. Line 241 Prelim.f"""
        return afiel_1 - awilt_1

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            pp_reg['afiel_1_path'], pp_reg['awilt_1_path']]],
        calc_wc, pp_reg['wc_path'], gdal.GDT_Float32, _TARGET_NODATA)

    def calc_eftext(peftxa, peftxb, sand):
        """Calculate effect of soil texture on microbial decomposition.

        Use an empirical regression to estimate the effect of soil
        sand content on the microbe decomposition rate. Line 359 Prelim.f

        Parameters:
            peftxa (numpy.ndarray): parameter, regression intercept
            peftxb (numpy.ndarray): parameter, regression slope
            sand (numpy.ndarray): proportion sand in soil

        Returns:
            eftext, coefficient that modifies microbe decomposition rate.
        """
        eftext = numpy.empty(sand.shape, dtype=numpy.float32)
        eftext[:] = _IC_NODATA
        valid_mask = (
            (peftxa != _TARGET_NODATA)
            & (peftxb != _TARGET_NODATA)
            & (sand != sand_nodata))
        eftext[valid_mask] = (
            peftxa[valid_mask] + (peftxb[valid_mask] * sand[valid_mask]))
        return eftext

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['peftxa'], param_val_dict['peftxb'], sand_path]],
        calc_eftext, pp_reg['eftext_path'], gdal.GDT_Float32, _IC_NODATA)

    def calc_p1co2_2(p1co2a_2, p1co2b_2, sand):
        """Calculate the fraction of carbon lost to CO2 from som1c_2.

        During decomposition from active organic soil carbon, a fraction
        of decomposing material is lost to CO2 as the soil respires.
        Line 366 Prelim.f

        Parameters:
            p1co2a_2 (numpy.ndarray): parameter, intercept of regression
                predicting loss to CO2 from active organic soil carbon
            p1co2b_2 (numpy.ndarray): parameter, slope of regression
                predicting loss to CO2 from active organic soil carbon
            sand (numpy.ndarray): proportion sand in soil

        Returns:
            p1co2_2, fraction of carbon that flows to CO2 from active
            organic soil carbon
        """
        p1co2_2 = numpy.empty(sand.shape, dtype=numpy.float32)
        p1co2_2[:] = _IC_NODATA
        valid_mask = (
            (p1co2a_2 != _TARGET_NODATA)
            & (p1co2b_2 != _TARGET_NODATA)
            & (sand != sand_nodata))
        p1co2_2[valid_mask] = (
            p1co2a_2[valid_mask] + (p1co2b_2[valid_mask] * sand[valid_mask]))
        return p1co2_2

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['p1co2a_2'], param_val_dict['p1co2b_2'], sand_path]],
        calc_p1co2_2, pp_reg['p1co2_2_path'], gdal.GDT_Float32, _IC_NODATA)

    def calc_fps1s3(ps1s3_1, ps1s3_2, clay):
        """Calculate effect of clay content on decomposition from som1c_2.

        Use an empirical regression to estimate the effect of clay content
        of soil on flow from soil organic matter with fast turnover to
        soil organic matter with slow turnover. Line 370 Prelim.f

        Parameters:
            ps1s3_1 (numpy.ndarray): parameter, regression intercept
            ps1s3_2 (numpy.ndarray): parameter, regression slope
            clay (numpy.ndarray): proportion clay in soil

        Returns:
            fps1s3, coefficient that modifies rate of decomposition
            from som1c_2
        """
        fps1s3 = numpy.empty(clay.shape, dtype=numpy.float32)
        fps1s3[:] = _IC_NODATA
        valid_mask = (
            (ps1s3_1 != _TARGET_NODATA)
            & (ps1s3_2 != _TARGET_NODATA)
            & (clay != clay_nodata))
        fps1s3[valid_mask] = (
            ps1s3_1[valid_mask] + (ps1s3_2[valid_mask] * clay[valid_mask]))
        return fps1s3

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['ps1s3_1'], param_val_dict['ps1s3_2'], clay_path]],
        calc_fps1s3, pp_reg['fps1s3_path'], gdal.GDT_Float32, _IC_NODATA)

    def calc_fps2s3(ps2s3_1, ps2s3_2, clay):
        """Calculate effect of clay content on decomposition from som2c_2.

        Use an empirical regression to estimate the effect of clay content
        of soil on flow from slow soil organic carbon to soil passive organic
        carbon. Line 371 Prelim.f

        Parameters:
            ps2s3_1 (numpy.ndarray): parameter, regression intercept
            ps2s3_2 (numpy.ndarray): parameter, regression slope
            clay (numpy.ndarray): proportion clay in soil

        Returns:
            fps2s3, coefficient that modifies rate of decomposition from
            som2c_2 to som3c
        """
        fps2s3 = numpy.empty(clay.shape, dtype=numpy.float32)
        fps2s3[:] = _IC_NODATA
        valid_mask = (
            (ps2s3_1 != _TARGET_NODATA)
            & (ps2s3_2 != _TARGET_NODATA)
            & (clay != clay_nodata))
        fps2s3[valid_mask] = (
            ps2s3_1[valid_mask] + (ps2s3_2[valid_mask] * clay[valid_mask]))
        return fps2s3

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['ps2s3_1'], param_val_dict['ps2s3_2'], clay_path]],
        calc_fps2s3, pp_reg['fps2s3_path'], gdal.GDT_Float32, _IC_NODATA)

    def calc_orglch(omlech_1, omlech_2, sand):
        """Calculate the effect of sand content on leaching from soil.

        Use an empirical regression to estimate the effect of sand content
        of soil on rate of organic leaching from soil when there is drainage
        of soil water from soil layer 1 to soil layer 2. Line 110 Predec.f

        Parameters:
            omlech_1 (numpy.ndarray): parameter, regression intercept
            omlech_2 (numpy.ndarray): parameter, regression slope
            sand (numpy.ndarray): proportion sand in soil

        Returns:
            orglch, the fraction of organic compounds leaching from soil
            with drainage from soil layer 1 to layer 2
        """
        orglch = numpy.empty(sand.shape, dtype=numpy.float32)
        orglch[:] = _IC_NODATA
        valid_mask = (
            (omlech_1 != _TARGET_NODATA)
            & (omlech_2 != _TARGET_NODATA)
            & (sand != sand_nodata))
        orglch[valid_mask] = (
            omlech_1[valid_mask] + (omlech_2[valid_mask] * sand[valid_mask]))
        return orglch

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['omlech_1'], param_val_dict['omlech_2'],
            sand_path]],
        calc_orglch, pp_reg['fps2s3_path'], gdal.GDT_Float32, _IC_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)

def _aboveground_ratio(anps, tca, pcemic_1, pcemic_2, pcemic_3, cemicb):
    """General function to calculate C/<iel> ratios of decomposing
    aboveground material. In general if these ratios are exceeeded, the
    material cannot decompose. Agdrat.f
    Inputs:
        anps = N or P in the donor material
        tca = total C in the donor material
        pcemic = fixed parameters
            pcemic_1 = maximum C/<iel> of new material
            pcemic_2 = minimum C/<iel> of new material
            pcemic_3 = minimum <iel> content of decomposing
                       material that gives minimum C/<iel> of new material
        cemicb = slope of the regression line for C/<iel>
    Returns:
        agdrat, the C/<iel> ratio of new material"""

    econt = np.where(tca > 0., anps / (tca * 2.5), 0.)
    agdrat = np.where(econt > pcemic_3, pcemic_2, pcemic_1 + econt * cemicb)
    return agdrat

def _structural_ratios(site_index_path, site_param_table, sv_reg, pp_reg):
    """Calculate maximum C/N and C/P ratios for decomposition of structural
    material (i.e., material containing lignin). These ratios do not change
    throughout the course of the simulation. Lines 31-77 Predec.f"""

    # temporary parameter rasters for these calculations
    temp_dir = tempfile.mkdtemp()
    param_val_dict = {}
    for iel in [1, 2]:
        for val in['pcemic1_2', 'pcemic1_1', 'pcemic1_3', 'pcemic2_2',
                   'pcemic2_1', 'pcemic2_3', 'rad1p_1', 'rad1p_2',
                   'rad1p_3', 'varat1_1', 'varat22_1']:
            target_path = os.path.join(temp_dir, '{}_{}.tif'.format(val, iel))
            param_val_dict['{}_{}'.format(val, iel)] = target_path
            site_to_val = dict(
                [(site_code, float(table['{}_{}'.format(val, iel)])) for
                (site_code, table) in site_param_table.items()])
            pygeoprocessing.reclassify_raster(
                (site_index_path, 1), site_to_val, target_path,
                gdal.GDT_Float32, _TARGET_NODATA)

    def calc_rnewas_som1(pcemic1_2, pcemic1_1, pcemic1_3,
                         struce_1, strucc_1):
        """Calculate maximum ratios for decomposition of aboveground
        structural material into SOM1."""
        valid_mask = strucc_1 >= 0.
        cemicb1 = np.empty(strucc_1.shape)
        cemicb1[:] = _TARGET_NODATA
        cemicb1[valid_mask] = (pcemic1_2[valid_mask] -
            pcemic1_1[valid_mask] / pcemic1_3[valid_mask])

        rnewas1 = _aboveground_ratio(struce_1, strucc_1, pcemic1_1,
            pcemic1_2, pcemic1_3, cemicb1)
        return rnewas1

    def calc_rnewas_som2(pcemic2_2, pcemic2_1, pcemic2_3,
                         struce_1, strucc_1,
                         rad1p_1, rad1p_2, rad1p_3,
                         pcemic1_2, rnewas1):
        """Calculate maximum ratios for decomposition of aboveground
        structural material into SOM2, including a fraction of the ratios
        entering SOM1."""
        valid_mask = strucc_1 >= 0.
        cemicb2 = np.empty(strucc_1.shape)
        cemicb2[:] = _TARGET_NODATA
        cemicb2[valid_mask] = (pcemic2_2[valid_mask] -
            pcemic2_1[valid_mask] / pcemic2_3[valid_mask])

        rnewas2 = _aboveground_ratio(struce_1, strucc_1, pcemic2_1,
            pcemic2_2, pcemic2_3, cemicb2)

        radds1 = np.empty(strucc_1.shape)
        radds1[:] = _TARGET_NODATA
        radds1[valid_mask] = (rad1p_1[valid_mask] + rad1p_2[valid_mask] *
            (rnewas1[valid_mask] - pcemic1_2[valid_mask]))
        rnewas2[valid_mask] = rnewas1[valid_mask] + radds1[valid_mask]
        rnewas2[valid_mask] = np.maximum(rnewas2[valid_mask],
            rad1p_3[valid_mask])
        return rnewas2

    for iel in [1, 2]:
        # calculate rnewas_iel_1 - aboveground material to SOM1
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                param_val_dict['pcemic1_2_{}'.format(iel)],
                param_val_dict['pcemic1_1_{}'.format(iel)],
                param_val_dict['pcemic1_3_{}'.format(iel)],
                sv_reg['struce_1_{}_path'.format(iel)],
                sv_reg['strucc_1_path']]],
            calc_rnewas_som1, pp_reg['rnewas_1_{}_path'.format(iel)],
            gdal.GDT_Float32, _TARGET_NODATA)
        # calculate rnewas_iel_2 - aboveground material to SOM2
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                param_val_dict['pcemic2_2_{}'.format(iel)],
                param_val_dict['pcemic2_1_{}'.format(iel)],
                param_val_dict['pcemic2_3_{}'.format(iel)],
                sv_reg['struce_1_{}_path'.format(iel)],
                sv_reg['strucc_1_path'],
                param_val_dict['rad1p_1_{}'.format(iel)],
                param_val_dict['rad1p_2_{}'.format(iel)],
                param_val_dict['rad1p_3_{}'.format(iel)],
                param_val_dict['pcemic1_2_{}'.format(iel)],
                pp_reg['rnewas_1_{}_path'.format(iel)]]],
            calc_rnewas_som2, pp_reg['rnewas_2_{}_path'.format(iel)],
            gdal.GDT_Float32, _TARGET_NODATA)
        # calculate rnewbs_iel_1 - belowground material to SOM1
        site_to_varat1_1 = dict([
            (site_code, float(table['varat1_1_{}'.format(iel)])) for
            (site_code, table) in site_param_table.items()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_varat1_1,
            pp_reg['rnewbs_{}_1_path'.format(iel)],
            gdal.GDT_Float32, _TARGET_NODATA)
        # calculate rnewbs_iel_2 - belowground material to SOM2
        # rnewbs(iel,2) = varat22(1,iel)
        site_to_varat22_1 = dict([
            (site_code, float(table['varat22_1_{}'.format(iel)])) for
            (site_code, table) in site_param_table.items()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_varat22_1,
            pp_reg['rnewbs_{}_2_path'.format(iel)],
            gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)

def _yearly_tasks(site_index_path, site_param_table, aligned_inputs,
                  month_index, year_reg):
    """Calculate annual precipitation and annual atmospheric N deposition.
    Century also calculates non-symbiotic soil N fixation once yearly, but
    here those were moved to monthly tasks.
    Century uses precipitation in the future 12 months (prcgrw) to predict
    root:shoot ratios, but here we use annual precipitation in 12 months
    including the current one instead if data for 12 future months are not
    available.
    Lines 79-82 Eachyr.f"""

    # annual precipitation = precip in 12 months following this one, or
    # in 12 months including months previous to this one if data for 12
    # future months are not available
    annual_precip_rasters = []
    for precip_month in xrange(month_index, month_index + 12):
        try:
            annual_precip_rasters.append(
                aligned_inputs['precip_%d' % precip_month])
        except KeyError:
            continue
    offset = 1
    while len(annual_precip_rasters) < 12:
        precip_month = month_index - offset
        try:
            annual_precip_rasters.append(
                    aligned_inputs['precip_%d' % precip_month])
        except KeyError:
            raise KeyError("Insufficient precipitation rasters were found")
        offset = offset + 1

    # sum them
    def raster_sum(*raster_list):
        return np.sum(raster_list, axis=0)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in annual_precip_rasters],
        raster_sum, year_reg['annual_precip_path'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate base N deposition
    # intermediate parameter rasters for this operation
    temp_dir = tempfile.mkdtemp()
    param_val_dict = {}
    for val in['epnfa_1', 'epnfa_2']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
            (site_code, table) in site_param_table.items()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_val, target_path,
            gdal.GDT_Float32, _TARGET_NODATA)

    def calc_base_N_dep(epnfa_1, epnfa_2, prcann):
        # baseNdep = max(epnfa_1 + epnfa_2 * MIN(prcann, 80.0), 0.)
        baseNdep = np.empty(prcann.shape)
        baseNdep[:] = 0.
        valid_mask = prcann >= 0.
        baseNdep[valid_mask] = epnfa_1[valid_mask] + (epnfa_2[valid_mask] *
            np.minimum(prcann[valid_mask], 80.))
        baseNdep[baseNdep < 0] = 0.
        return baseNdep

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [param_val_dict['epnfa_1'],
            param_val_dict['epnfa_2'], year_reg['annual_precip_path']]],
        calc_base_N_dep, year_reg['baseNdep_path'], gdal.GDT_Float32,
        _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)

def _shortwave_radiation(template_raster, month, shwave_path):
    """Calculate shortwave radiation outside the atmosphere, using the input
    `template_raster` to calculate the latitude of each pixel. shwave.f"""

    def shwave(month):
        def _shwave(latitude):
            """Inputs:
            latitude - latitude of current site in degrees
            month - current month

            Output:
            shwave - short wave solar radiation outside the atmosphere

            Local variables:
            ahou            - ?
            declin          - declination (radians)
            jday()          - Julian day for middle of current month
            par1, par2      - parameters in computation of ahou
            rlatitude       - latitude of the site (in radians)
            solrad          - solar radiation (ly/day)
            transcof        - transmission coefficient"""

            # Julian date in middle of each month of the year
            jday_list = [
                16, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350]
            jday = jday_list[month]

            transcof = 0.8
            valid_mask = latitude >= -180.

            rlatitude = numpy.empty(latitude.shape)
            rlatitude[:] = _IC_NODATA

            # Convert latitude from degrees to radians
            rlatitude[valid_mask] = latitude[valid_mask] * (numpy.pi / 180.0)

            # short wave solar radiation on a clear day
            declin = 0.401426 * numpy.sin(6.283185 * (jday - 77.0) / 365.0)

            temp = 1.0 - (-numpy.tan(rlatitude) * numpy.tan(declin))**2
            temp = numpy.where(temp < 0., 0., temp)

            par1 = numpy.sqrt(temp)
            par2 = (-numpy.tan(rlatitude) * numpy.tan(declin))

            ahou = numpy.arctan2(par1, par2)
            ahou = numpy.where(ahou < 0., 0., ahou)

            solrad = 917.0 * transcof * (ahou * numpy.sin(rlatitude) *
                numpy.sin(declin) + numpy.cos(rlatitude) * numpy.cos(declin) *
                numpy.sin(ahou))

            # short wave radiation outside the atmosphere
            shwave = numpy.empty(latitude.shape)
            shwave[:] = _TARGET_NODATA
            shwave[valid_mask] = solrad[valid_mask] / transcof
            return shwave
        return _shwave

    # TODO if we allow projected inputs in the future, must reproject the
    # template raster here to ensure we collect geographic coordinates
    # calculate an intermediate input, latitude at each pixel center
    temp_dir = tempfile.mkdtemp()
    latitude_raster_path = os.path.join(temp_dir, 'latitude.tif')
    pygeoprocessing.new_raster_from_base(
        template_raster, latitude_raster_path, gdal.GDT_Float32,
        [_IC_NODATA])
    latitude_raster = gdal.OpenEx(latitude_raster_path, gdal.GA_Update)
    target_band = latitude_raster.GetRasterBand(1)
    base_raster_info = pygeoprocessing.get_raster_info(template_raster)
    geotransform = base_raster_info['geotransform']
    for offset_map, raster_block in pygeoprocessing.iterblocks(
            template_raster):
        n_y_block = raster_block.shape[0]
        n_x_block = raster_block.shape[1]

        # offset by .5 so we're in the center of the pixel
        xoff = offset_map['xoff'] + 0.5
        yoff = offset_map['yoff'] + 0.5

        # calculate the projected x and y coordinate bounds for the block
        x_range = numpy.linspace(
            geotransform[0] + geotransform[1] * xoff,
            geotransform[0] + geotransform[1] * (xoff + n_x_block - 1),
            n_x_block)
        y_range = numpy.linspace(
            geotransform[3] + geotransform[5] * yoff,
            geotransform[3] + geotransform[5] * (yoff + n_y_block - 1),
            n_y_block)

        # we'll use this to avoid generating any nodata points
        valid_mask = raster_block != base_raster_info['nodata']

        # these indexes correspond to projected coordinates
        # y_vector is what we want, an array of latitude coordinates
        x_vector, y_vector = numpy.meshgrid(x_range, y_range)

        target_band.WriteArray(y_vector, xoff=offset_map['xoff'],
            yoff=offset_map['yoff'])

    # Making sure the band and dataset is flushed and not in memory
    target_band.FlushCache()
    target_band.FlushCache()
    target_band = None
    gdal.Dataset.__swig_destroy__(latitude_raster)
    latitude_raster = None

    pygeoprocessing.raster_calculator(
        [(latitude_raster_path, 1)],
        shwave(month), shwave_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)

def _reference_evaporation(max_temp, min_temp, shwave, fwloss_4):
    """Calculate reference evapotranspiration with the FAO Penman-Monteith
    equation (http://www.fao.org/docrep/X0490E/x0490e08.htm),
    modified by the parameter fwloss(4). Pevap.f"""

    const1 = 0.0023
    const2 = 17.8
    langleys2watts = 54.0

    valid_mask = (fwloss_4 > 0) & (shwave >= 0)
    trange = np.empty(fwloss_4.shape)
    trange[:] = _TARGET_NODATA
    trange[valid_mask] = max_temp[valid_mask] - min_temp[valid_mask]
    tmean = np.empty(fwloss_4.shape)
    tmean[:] = _IC_NODATA
    tmean[valid_mask] = (max_temp[valid_mask] + min_temp[valid_mask]) / 2.0

    # daily reference evapotranspiration
    daypet = np.empty(fwloss_4.shape)
    daypet[:] = _TARGET_NODATA
    in1 = const1 * (tmean[valid_mask] + const2)
    in2 = np.sqrt(trange[valid_mask])
    in3 = (shwave[valid_mask] / langleys2watts)
    daypet[valid_mask] = (const1 * (tmean[valid_mask] + const2) *
        np.sqrt(trange[valid_mask]) * (shwave[valid_mask] / langleys2watts))

    # monthly reference evapotranspiration, from mm to cm,
    # bounded to be at least 0.5
    monpet = np.where(((daypet * 30.) / 10.) > 0.5,
        ((daypet * 30.) / 10.), 0.5)

    pevap = np.empty(fwloss_4.shape)
    pevap[:] = _TARGET_NODATA
    pevap[valid_mask] = monpet[valid_mask] * fwloss_4[valid_mask]
    return pevap

def _potential_production(aligned_inputs, site_param_table, current_month,
                          pft_id_set, veg_trait_table, sv_reg):
    """Calculate total potential production and root:shoot ratios. Potcrp.f"""
    # if growth does not occur this month for all PFTs,
    # skip the rest of the function
    do_PFT = []
    for pft_index in pft_id_set:
        if str(current_month) in veg_trait_table[pft_index]['growth_months']:
            do_PFT.append(pft_index)
    if not do_PFT:
        return

    def multiply_positive_rasters(raster1, raster2):
        masked_r1 = np.where(raster1 < 0., 0., raster1)
        masked_r2 = np.where(raster2 < 0., 0., raster2)
        return masked_r1 * masked_r2

    def sum_positive_rasters(*raster_list):
        masked_list = [np.where(r < 0., 0., r) for r in raster_list]
        return np.sum(masked_list, axis=0)

    def calc_ctemp(aglivc, pmxbio, maxtmp, pmxtmp, mintmp, pmntmp):
        """Calculate soil temperature relative to its effect on growth.
        Lines 69-84 Potcrp.f"""
        bio = np.empty(aglivc.shape)
        bio[:] = _IC_NODATA
        valid_mask = pmxbio >= 0
        bio[valid_mask] = aglivc[valid_mask] * 2.5
        bio = np.where(bio > pmxbio, pmxbio, bio)
        bio[pmxbio < 0] = _IC_NODATA

        # Maximum temperature
        tmxs = np.empty(aglivc.shape)
        tmxs[:] = _IC_NODATA
        tmxs[valid_mask] = maxtmp[valid_mask] + ((25.4/(1. + 18. *
            np.exp(-0.20 * maxtmp[valid_mask]))) *
            (np.exp(pmxtmp[valid_mask] * bio[valid_mask]) - 0.13))

        # Minimum temperature
        tmns = np.empty(aglivc.shape)
        tmns[:] = _IC_NODATA
        tmns[valid_mask] = mintmp[valid_mask] + (pmntmp[valid_mask] *
            bio[valid_mask] - 1.78)

        # Average temperature
        ctemp = np.empty(aglivc.shape)
        ctemp[:] = _IC_NODATA
        ctemp[valid_mask] = (tmxs[valid_mask] + tmns[valid_mask])/2.
        return ctemp

    # temporary intermediate rasters for these calculations
    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {}
    # site-level temporary calculated values
    for val in ['sum_aglivc', 'sum_stdedc', 'ctemp', 'shwave', 'pevap']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))
    # PFT-level temporary calculated values
    for pft_i in pft_id_set:
        for val in ['aglivc_weighted', 'stdedc_weighted']:
            temp_val_dict['{}_{}'.format(val, pft_i)] = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))

    # temporary parameter rasters for these calculations
    param_val_dict = {}
    # site-level parameters
    for val in ['pmxbio', 'pmxtmp', 'pmntmp', 'fwloss_4']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
            (site_code, table) in site_param_table.items()])
        pygeoprocessing.reclassify_raster(
            (aligned_inputs['site_index'], 1), site_to_val, target_path,
            gdal.GDT_Float32, _TARGET_NODATA)
    # PFT-level parameters

    # calculate intermediate quantities that do not differ between PFTs:
    # sum of aglivc (standing live biomass) and stdedc (standing dead biomass)
    # across PFTs, weighted by % cover of each PFT
    for sv in ['aglivc', 'stdedc']:
        weighted_path_list = []
        for pft_index in pft_id_set:
            target_path = temp_val_dict['{}_weighted_{}'.format(sv, pft_index)]
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in
                    sv_reg['{}_{}_path'.format(sv, pft_index)],
                    aligned_inputs['pft_{}'.format(pft_index)]],
                multiply_positive_rasters, target_path,
                gdal.GDT_Float32, _TARGET_NODATA)
            weighted_path_list.append(target_path)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in weighted_path_list],
            sum_positive_rasters, temp_val_dict['sum_{}'.format(sv)],
            gdal.GDT_Float32, _TARGET_NODATA)

    # ctemp, soil temperature relative to impacts on growth
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [temp_val_dict['sum_aglivc'],
            param_val_dict['pmxbio'],
            aligned_inputs['max_temp_{}'.format(current_month)],
            param_val_dict['pmxtmp'],
            aligned_inputs['min_temp_{}'.format(current_month)],
            param_val_dict['pmntmp']]],
        calc_ctemp, temp_val_dict['ctemp'], gdal.GDT_Float32, _IC_NODATA)

    # shwave, shortwave radiation outside the atmosphere
    _shortwave_radiation(aligned_inputs['site_index'], current_month,
        temp_val_dict['shwave'])

    # pet, reference evapotranspiration modified by fwloss parameter
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            aligned_inputs['max_temp_{}'.format(current_month)],
            aligned_inputs['min_temp_{}'.format(current_month)],
            temp_val_dict['shwave'],
            param_val_dict['fwloss_4']]],
        _reference_evaporation, temp_val_dict['pevap'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # h2ogef_prior, line 59 Potcrp.f

    # calculate quantities that differ between PFTs
    # for pft_i in do_PFT:
        # potprd, the limiting effect of temperature on growth
        # potprd = f(temp_val_dict['ctemp'])
        # h2ogef(1), the limiting effect of soil water availability on growth
        # h2ogef(1) = f(pet) # does this need to take into account other PFTs?
        # biof, the limiting effect of obstruction by standing dead on growth
        # biof = f(sum of STDEDC, STRUCC_1) sv_reg['strucc_1_path']
        # sdlng, the limiting effect of shading by live material on growth
        # sdlng = f(sum of AGLIVC)
        # total potential production for each PFT
        # tgprod = prdx(1) * shwave * potprd * h2ogef(1) * biof * sdlng

    # clean up temporary files
    shutil.rmtree(temp_dir)
