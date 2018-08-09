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
    'minerl_1_1_path': 'minerl_1_1.tif',
    'minerl_2_1_path': 'minerl_2_1.tif',
    'minerl_3_1_path': 'minerl_3_1.tif',
    'minerl_4_1_path': 'minerl_4_1.tif',
    'minerl_5_1_path': 'minerl_5_1.tif',
    'minerl_6_1_path': 'minerl_6_1.tif',
    'minerl_7_1_path': 'minerl_7_1.tif',
    'minerl_8_1_path': 'minerl_8_1.tif',
    'minerl_9_1_path': 'minerl_9_1.tif',
    'minerl_10_1_path': 'minerl_10_1.tif',
    'minerl_1_2_path': 'minerl_1_2.tif',
    'minerl_2_2_path': 'minerl_2_2.tif',
    'minerl_3_2_path': 'minerl_3_2.tif',
    'minerl_4_2_path': 'minerl_4_2.tif',
    'minerl_5_2_path': 'minerl_5_2.tif',
    'minerl_6_2_path': 'minerl_6_2.tif',
    'minerl_7_2_path': 'minerl_7_2.tif',
    'minerl_8_2_path': 'minerl_8_2.tif',
    'minerl_9_2_path': 'minerl_9_2.tif',
    'minerl_10_2_path': 'minerl_10_2.tif',
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
    'stdede_1', 'aglive_2', 'bglive_2', 'stdede_2', 'avh2o_1',
    'crpstg_1', 'crpstg_2',
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

# intermediate values for each plant functional type that are shared
# between submodels, but do not need to be saved as output
_PFT_INTERMEDIATE_VALUES = [
    'h2ogef_1', 'tgprod_pot_prod',
    'cercrp_min_above_1', 'cercrp_min_above_2',
    'cercrp_max_above_1', 'cercrp_max_above_2',
    'cercrp_min_below_1', 'cercrp_min_below_2',
    'cercrp_max_below_1', 'cercrp_max_below_2',
    'tgprod', 'rtsh']

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
            "Couldn't find the following precipitation paths given the " +
            "pattern: %s\n\t" % args['monthly_precip_path_pattern'] +
            "\n\t".join(missing_precip_path_list))
    # the model requires 12 months of precipitation data to calculate
    # atmospheric N deposition and potential production from annual precip
    n_precip_months = int(args['n_months'])
    if n_precip_months < 12:
        m_index = int(args['n_months'])
        while m_index <= 12:
            month_i = (starting_month + m_index - 1) % 12 + 1
            year = starting_year + (starting_month + m_index - 1) // 12
            precip_path = args['monthly_precip_path_pattern'].replace(
                '<year>', str(year)).replace('<month>', '%.2d' % month_i)
            base_align_raster_path_id_map['precip_%d' % m_index] = precip_path
            precip_path_list.append(precip_path)
            if os.path.exists(precip_path):
                n_precip_months = n_precip_months + 1
            m_index = m_index + 1
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
                "Couldn't find the following temperature raster paths" +
                " given the pattern: %s\n\t" % args[
                    '%s_temp_path_pattern' % substring] +
                "\n\t".join(missing_temperature_path_list))

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
            "Couldn't find parameter values for the following site " +
            "indices: %s\n\t" + ", ".join(missing_site_index_list))

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
            "Couldn't find trait values for the following plant functional " +
            "types: %s\n\t" + ", ".join(missing_pft_trait_list))
    frtcindx_set = set([
        veg_trait_table[k]['frtcindx'] for k in veg_trait_table.iterkeys()])
    if frtcindx_set.difference(set([0, 1])):
        raise ValueError("frtcindx parameter contains invalid values")

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
            "Couldn't find trait values for the following animal " +
            "ids: %s\n\t" + ", ".join(missing_animal_trait_list))

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
                "Couldn't find the following required initial values: " +
                "\n\t".join(missing_initial_values))

        # align and resample initialization rasters with inputs
        sv_dir = os.path.join(args['workspace_dir'], 'state_variables_m-1')
        aligned_initial_path_map = dict(
            [(key, os.path.join(sv_dir, os.path.basename(path)))
                for key, path in resample_initial_path_map.iteritems()])
        initial_path_list = (
            [resample_initial_path_map[k] for k in
                sorted(resample_initial_path_map.iterkeys())] +
            source_input_path_list)
        aligned_initial_path_list = (
            [aligned_initial_path_map[k] for k in
                sorted(aligned_initial_path_map.iterkeys())] +
            aligned_input_path_list)
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

    # calculate required ratios for decomposition of structural material
    LOGGER.info("Calculating required ratios for structural decomposition")
    _structural_ratios(
        aligned_inputs['site_index'], site_param_table, sv_reg, pp_reg)

    # make yearly directory for values that are updated every twelve months
    year_dir = tempfile.mkdtemp()
    year_reg = dict(
        [(key, os.path.join(year_dir, path)) for key, path in
            _YEARLY_FILES.iteritems()])

    # make monthly directory for monthly intermediate parameters that are
    # shared between submodels, but do not need to be saved as output
    month_temp_dir = tempfile.mkdtemp()
    month_reg = {}
    for pft_index in pft_id_set:
        for val in _PFT_INTERMEDIATE_VALUES:
            month_reg['{}_{}'.format(
                val, pft_index)] = os.path.join(
                month_temp_dir, '{}_{}.tif'.format(val, pft_index))

    # Main simulation loop
    # for each step in the simulation
    for month_index in xrange(n_months):
        if (month_index % 12) == 0:
            # Update yearly quantities
            _yearly_tasks(
                aligned_inputs['site_index'], site_param_table,
                aligned_inputs, month_index, year_reg)

        current_month = (starting_month + month_index - 1) % 12 + 1
        year = starting_year + (starting_month + month_index - 1) // 12

        # make new folders for state variables during this step
        sv_dir = os.path.join(
            args['workspace_dir'], 'state_variables_m%d' % month_index)
        utils.make_directories([sv_dir])

        # track state variables from previous step
        # prev_sv_reg = sv_reg
        # sv_reg = utils.build_file_registry(
        #     [(_SITE_STATE_VARIABLE_FILES, sv_dir),
        #         (pft_sv_dict, sv_dir)], file_suffix)

        # update state variables from previous month
        LOGGER.info(
            "Main simulation loop: month %d of %d" % (
                month_index, n_months))

        _potential_production(
            aligned_inputs, site_param_table, current_month, month_index,
            pft_id_set, veg_trait_table, sv_reg, pp_reg, month_reg)

        _root_shoot_ratio(
            aligned_inputs, site_param_table, current_month, pft_id_set,
            veg_trait_table, sv_reg, year_reg, month_reg)


def raster_sum(
        raster_list, input_nodata, target_path, target_nodata,
        nodata_remove=False):
    """Calculate the sum per pixel across rasters in a list.

    Sum the rasters in `raster_list` element-wise, allowing nodata values
    in the rasters to propagate to the result or treating nodata as zero.

    Parameters:
        raster_list (list): list of paths to rasters to sum
        input_nodata (float or int): nodata value in the input rasters
        target_path (string): path to location to store the result
        target_nodata (float or int): nodata value for the result raster
        nodata_remove (bool): if true, treat nodata values in input
            rasters as zero. If false, the sum in a pixel where any input
            raster is nodata is nodata.

    Modifies:
        the raster indicated by `target_path`

    Returns:
        None
    """
    def raster_sum_op(*raster_list):
        """Add the rasters in raster_list without removing nodata values."""
        invalid_mask = numpy.any(
            numpy.array(raster_list) == input_nodata, axis=0)
        sum_of_rasters = numpy.sum(raster_list, axis=0)
        sum_of_rasters[invalid_mask] = target_nodata
        return sum_of_rasters

    def raster_sum_op_nodata_remove(*raster_list):
        """Add the rasters in raster_list, treating nodata as zero."""
        for r in raster_list:
            numpy.place(r, r == input_nodata, [0])
        sum_of_rasters = numpy.sum(raster_list, axis=0)
        return sum_of_rasters

    if nodata_remove:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in raster_list], raster_sum_op_nodata_remove,
            target_path, gdal.GDT_Float32, target_nodata)

    else:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in raster_list], raster_sum_op,
            target_path, gdal.GDT_Float32, target_nodata)


def _calc_ompc(
        som1c_2_path, som2c_2_path, som3c_path, bulkd_path, edepth_path,
        ompc_path):
    """Estimate total soil organic matter.

    Total soil organic matter is the sum of soil carbon across
    slow, active, and passive compartments, weighted by bulk
    density and total modeled soil depth. Lines 220-222, Prelim.f

    Parameters:
        som1c_2_path (string): path to active organic soil carbon raster
        som2c_2_path (string): path to slow organic soil carbon raster
        som3c_path (string): path to passive organic soil carbon raster
        bulkd_path (string): path to bulk density of soil raster
        edepth (string): path to depth of soil raster
        ompc_path (string): path to result, total soil organic matter

    Modifies:
        the raster indicated by `ompc_path`

    Returns:
        None
    """
    def ompc_op(som1c_2, som2c_2, som3c, bulkd, edepth):
        """Estimate total soil organic matter.

        Total soil organic matter is the sum of soil carbon across
        slow, active, and passive compartments, weighted by bulk
        density and total modeled soil depth. Lines 220-222, Prelim.f

        Parameters:
            som1c_2_path (string): state variable, active organic soil carbon
            som2c_2_path (string): state variable, slow organic soil carbon
            som3c_path (string): state variable, passive organic soil carbon
            bulkd_path (string): input, bulk density of soil
            edepth_path (string): parameter, depth of soil for this
                calculation

        Returns:
            ompc, total soil organic matter weighted by bulk
            density.
        """
        ompc = numpy.empty(som1c_2.shape, dtype=numpy.float32)
        ompc[:] = _TARGET_NODATA
        valid_mask = (
            (som1c_2 != som1c_2_nodata) &
            (som2c_2 != som2c_2_nodata) &
            (som3c != som3c_nodata) &
            (bulkd != bulkd_nodata) &
            (edepth != _IC_NODATA))
        ompc[valid_mask] = (
            (som1c_2[valid_mask] + som2c_2[valid_mask] +
                som3c[valid_mask]) * 1.724 /
            (10000. * bulkd[valid_mask] * edepth[valid_mask]))
        return ompc

    som1c_2_nodata = pygeoprocessing.get_raster_info(som1c_2_path)['nodata'][0]
    som2c_2_nodata = pygeoprocessing.get_raster_info(som2c_2_path)['nodata'][0]
    som3c_nodata = pygeoprocessing.get_raster_info(som3c_path)['nodata'][0]
    bulkd_nodata = pygeoprocessing.get_raster_info(bulkd_path)['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            som1c_2_path, som2c_2_path, som3c_path,
            bulkd_path, edepth_path]],
        ompc_op, ompc_path, gdal.GDT_Float32, _TARGET_NODATA)


def _calc_afiel(
        sand_path, silt_path, clay_path, ompc_path, bulkd_path, afiel_path):
    """Calculate field capacity for one soil layer.

    Parameters:
        sand_path (string): path to proportion sand in soil raster
        silt_path (string): path to proportion silt in soil raster
        clay_path (string): path to proportion clay in soil raster
        ompc_path (string): path to estimated total soil organic matter raster
        bulkd_path (string): path to bulk density of soil raster
        afiel_path (string): path to result raster, field capacity for this
            soil layer

    Modifies:
        the raster indicated by `afiel_path`

    Returns:
        None
    """
    def afiel_op(sand, silt, clay, ompc, bulkd):
        """Calculate field capacity for one soil layer.

        Field capacity, maximum soil moisture retention capacity,
        from Gupta and Larson 1979, 'Estimating soil and water
        retention characteristics from particle size distribution,
        organic matter percent and bulk density'. Water Resources
        Research 15:1633.

        Parameters:
            sand_path (string): input, proportion sand in soil
            silt_path (string): input, proportion silt in soil
            clay_path (string): input, proportion clay in soil
            ompc_path (string): derived, estimated total soil organic matter
            bulkd_path (string): input, bulk density of soil

        Returns:
            afiel, field capacity for this soil layer
        """
        afiel = numpy.empty(sand.shape, dtype=numpy.float32)
        afiel[:] = _TARGET_NODATA
        valid_mask = (
            (sand != sand_nodata) &
            (silt != silt_nodata) &
            (clay != clay_nodata) &
            (ompc != _TARGET_NODATA) &
            (bulkd != bulkd_nodata))
        afiel[valid_mask] = (
            0.3075 * sand[valid_mask] + 0.5886 * silt[valid_mask] +
            0.8039 * clay[valid_mask] + 2.208E-03 * ompc[valid_mask] +
            -0.1434 * bulkd[valid_mask])
        return afiel

    sand_nodata = pygeoprocessing.get_raster_info(sand_path)['nodata'][0]
    silt_nodata = pygeoprocessing.get_raster_info(silt_path)['nodata'][0]
    clay_nodata = pygeoprocessing.get_raster_info(clay_path)['nodata'][0]
    bulkd_nodata = pygeoprocessing.get_raster_info(bulkd_path)['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            sand_path, silt_path, clay_path, ompc_path, bulkd_path]],
        afiel_op, afiel_path, gdal.GDT_Float32, _TARGET_NODATA)


def _calc_awilt(
        sand_path, silt_path, clay_path, ompc_path, bulkd_path, awilt_path):
    """Calculate wilting point for one soil layer.

    Wilting point, minimum soil water required by plants before
    wilting, from Gupta and Larson 1979, 'Estimating soil and
    water retention characteristics from particle size distribution,
    organic matter percent and bulk density'. Water Resources
    Research 15:1633.

    Parameters:
        sand_path (string): path to proportion sand in soil raster
        silt_path (string): path to proportion silt in soil raster
        clay_path (string): path to proportion clay in soil raster
        ompc_path (string): path to estimated total soil organic matter raster
        bulkd_path (string): path to bulk density of soil raster
        awilt_path (string): path to result raster, wilting point for this
            soil layer

    Modifies:
        the raster indicated by `awilt_path`

    Returns:
        None
    """
    def awilt_op(sand, silt, clay, ompc, bulkd):
        """Calculate wilting point for one soil layer.

        Wilting point, minimum soil water required by plants before
        wilting, from Gupta and Larson 1979, 'Estimating soil and
        water retention characteristics from particle size distribution,
        organic matter percent and bulk density'. Water Resources
        Research 15:1633.

        Parameters:
            sand_path (string): input, proportion sand in soil
            silt_path (string): input, proportion silt in soil
            clay_path (string): input, proportion clay in soil
            ompc_path (string): derived, estimated total soil organic matter
            bulkd_path (string): input, bulk density of soil

        Returns:
            awilt, wilting point for this soil layer
        """
        awilt = numpy.empty(sand.shape, dtype=numpy.float32)
        awilt[:] = _TARGET_NODATA
        valid_mask = (
            (sand != sand_nodata) &
            (silt != silt_nodata) &
            (clay != clay_nodata) &
            (ompc != _TARGET_NODATA) &
            (bulkd != bulkd_nodata))
        awilt[valid_mask] = (
            -0.0059 * sand[valid_mask] + 0.1142 * silt[valid_mask] +
            0.5766 * clay[valid_mask] + 2.228E-03 * ompc[valid_mask] +
            0.02671 * bulkd[valid_mask])
        return awilt

    sand_nodata = pygeoprocessing.get_raster_info(sand_path)['nodata'][0]
    silt_nodata = pygeoprocessing.get_raster_info(silt_path)['nodata'][0]
    clay_nodata = pygeoprocessing.get_raster_info(clay_path)['nodata'][0]
    bulkd_nodata = pygeoprocessing.get_raster_info(bulkd_path)['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            sand_path, silt_path, clay_path, ompc_path, bulkd_path]],
        awilt_op, awilt_path, gdal.GDT_Float32, _TARGET_NODATA)


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
            the simulation

    Modifies the rasters pp_reg['afiel_<layer>'] and pp_reg['awilt_<layer>']
        for all soil layers.

    Returns:
        None
    """
    def decrement_ompc(ompc_orig_path, ompc_dec_path):
        """Decrease estimated organic matter to 85% of its value.

        In each subsequent soil layer, estimated organic matter is decreased
        by 15%, to 85% of its previous value.

        Parameters:
            ompc_orig_path (string): path to estimated soil organic matter
                raster
            ompc_dec_path (string): path to result raster, estimated soil
                organic matter decreased to 85% of its previous value

        Modifies:
            the raster indicated by `ompc_dec_path`

        Returns:
            None
        """
        def decrement_op(ompc_orig):
            """Reduce organic matter to 85% of its previous value."""
            ompc_dec = numpy.empty(ompc_orig.shape, dtype=numpy.float32)
            ompc_dec[:] = _TARGET_NODATA
            valid_mask = (ompc_orig != _TARGET_NODATA)
            ompc_dec[valid_mask] = ompc_orig[valid_mask] * 0.85
            return ompc_dec

        pygeoprocessing.raster_calculator(
            [(ompc_orig_path, 1)], decrement_op, ompc_dec_path,
            gdal.GDT_Float32, _TARGET_NODATA)

    # temporary intermediate rasters for calculating field capacity and
    # wilting point
    temp_dir = tempfile.mkdtemp()
    edepth_path = os.path.join(temp_dir, 'edepth.tif')
    ompc_path = os.path.join(temp_dir, 'ompc.tif')

    site_to_edepth = dict(
        [(site_code, float(table['edepth'])) for
         (site_code, table) in site_param_table.iteritems()])

    pygeoprocessing.reclassify_raster(
        (site_index_path, 1), site_to_edepth, edepth_path, gdal.GDT_Float32,
        _IC_NODATA)

    # estimate total soil organic matter
    _calc_ompc(
        som1c_2_path, som2c_2_path, som3c_path, bulk_d_path, edepth_path,
        ompc_path)

    # calculate field capacity and wilting point for each soil layer,
    # decreasing organic matter content by 85% with each layer
    for lyr in xrange(1, 10):
        afiel_path = pp_reg['afiel_{}_path'.format(lyr)]
        awilt_path = pp_reg['awilt_{}_path'.format(lyr)]
        _calc_afiel(
            sand_path, silt_path, clay_path, ompc_path, bulk_d_path,
            afiel_path)
        _calc_awilt(
            sand_path, silt_path, clay_path, ompc_path, bulk_d_path,
            awilt_path)
        ompc_dec_path = os.path.join(temp_dir, 'ompc{}.tif'.format(lyr))
        decrement_ompc(ompc_path, ompc_dec_path)
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
        pp_reg['wc_path']
        pp_reg['eftext_path']
        pp_reg['p1co2_2_path']
        pp_reg['fps1s3_path']
        pp_reg['fps2s3_path']
        pp_reg['orglch_path']

    Returns:
        None
    """
    sand_nodata = pygeoprocessing.get_raster_info(sand_path)['nodata'][0]
    clay_nodata = pygeoprocessing.get_raster_info(clay_path)['nodata'][0]

    # temporary intermediate rasters for persistent parameters calculation
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
            _IC_NODATA)

    def calc_wc(afiel_1, awilt_1):
        """Calculate water content of soil layer 1."""
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
            sand (numpy.ndarray): input, proportion sand in soil

        Returns:
            eftext, coefficient that modifies microbe decomposition rate.
        """
        eftext = numpy.empty(sand.shape, dtype=numpy.float32)
        eftext[:] = _IC_NODATA
        valid_mask = (
            (peftxa != _IC_NODATA) &
            (peftxb != _IC_NODATA) &
            (sand != sand_nodata))
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
            sand (numpy.ndarray): input, proportion sand in soil

        Returns:
            p1co2_2, fraction of carbon that flows to CO2 from active
            organic soil carbon
        """
        p1co2_2 = numpy.empty(sand.shape, dtype=numpy.float32)
        p1co2_2[:] = _IC_NODATA
        valid_mask = (
            (p1co2a_2 != _IC_NODATA) &
            (p1co2b_2 != _IC_NODATA) &
            (sand != sand_nodata))
        p1co2_2[valid_mask] = (
            p1co2a_2[valid_mask] + (p1co2b_2[valid_mask] * sand[valid_mask]))
        return p1co2_2

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['p1co2a_2'],
            param_val_dict['p1co2b_2'], sand_path]],
        calc_p1co2_2, pp_reg['p1co2_2_path'], gdal.GDT_Float32, _IC_NODATA)

    def calc_fps1s3(ps1s3_1, ps1s3_2, clay):
        """Calculate effect of clay content on decomposition from som1c_2.

        Use an empirical regression to estimate the effect of clay content
        of soil on flow from soil organic matter with fast turnover to
        soil organic matter with slow turnover. Line 370 Prelim.f

        Parameters:
            ps1s3_1 (numpy.ndarray): parameter, regression intercept
            ps1s3_2 (numpy.ndarray): parameter, regression slope
            clay (numpy.ndarray): input, proportion clay in soil

        Returns:
            fps1s3, coefficient that modifies rate of decomposition
            from som1c_2
        """
        fps1s3 = numpy.empty(clay.shape, dtype=numpy.float32)
        fps1s3[:] = _IC_NODATA
        valid_mask = (
            (ps1s3_1 != _IC_NODATA) &
            (ps1s3_2 != _IC_NODATA) &
            (clay != clay_nodata))
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
            clay (numpy.ndarray): input, proportion clay in soil

        Returns:
            fps2s3, coefficient that modifies rate of decomposition from
            som2c_2 to som3c
        """
        fps2s3 = numpy.empty(clay.shape, dtype=numpy.float32)
        fps2s3[:] = _IC_NODATA
        valid_mask = (
            (ps2s3_1 != _IC_NODATA) &
            (ps2s3_2 != _IC_NODATA) &
            (clay != clay_nodata))
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
            sand (numpy.ndarray): input, proportion sand in soil

        Returns:
            orglch, the fraction of organic compounds leaching from soil
            with drainage from soil layer 1 to layer 2
        """
        orglch = numpy.empty(sand.shape, dtype=numpy.float32)
        orglch[:] = _IC_NODATA
        valid_mask = (
            (omlech_1 != _IC_NODATA) &
            (omlech_2 != _IC_NODATA) &
            (sand != sand_nodata))
        orglch[valid_mask] = (
            omlech_1[valid_mask] + (omlech_2[valid_mask] * sand[valid_mask]))
        return orglch

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['omlech_1'], param_val_dict['omlech_2'],
            sand_path]],
        calc_orglch, pp_reg['orglch_path'], gdal.GDT_Float32, _IC_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _aboveground_ratio(anps, tca, pcemic_1, pcemic_2, pcemic_3, cemicb):
    """Calculate C/<iel> ratios of decomposing aboveground material.

    This ratio is used to test whether there is sufficient <iel> (N or P)
    in aboveground material for the material to decompose. Agdrat.f

    Parameters:
        anps (numpy.ndarray): state variable, N or P in the donor material
        tca (numpy.ndarray): state variable, total C in the donor material
        pcemic_1 (numpy.ndarray): parameter, maximum C/<iel> of new material
        pcemic_2 (numpy.ndarray): parameter, minimum C/<iel> of new material
        pcemic_3 (numpy.ndarray): parameter, minimum <iel> content of
            decomposing material that gives minimum C/<iel> of new material
        cemicb (numpy.ndarray): parameter, slope of the regression line for
            C/<iel>

    Returns:
        agdrat, the C/<iel> ratio of new material
    """
    valid_mask = (
        (anps != _TARGET_NODATA) &
        (tca != _TARGET_NODATA) &
        (pcemic_1 != _IC_NODATA) &
        (pcemic_2 != _IC_NODATA) &
        (pcemic_3 != _IC_NODATA) &
        (cemicb != _IC_NODATA))
    econt = numpy.empty(anps.shape, dtype=numpy.float32)
    econt[:] = _TARGET_NODATA
    econt[valid_mask] = numpy.where(
        tca[valid_mask] > 0., anps[valid_mask] / (tca[valid_mask] * 2.5), 0.)

    agdrat = numpy.empty(anps.shape, dtype=numpy.float32)
    agdrat[:] = _TARGET_NODATA
    agdrat[valid_mask] = numpy.where(
        econt[valid_mask] > pcemic_3[valid_mask],
        pcemic_2[valid_mask],
        pcemic_1[valid_mask] + econt[valid_mask] * cemicb[valid_mask])
    return agdrat


def _structural_ratios(site_index_path, site_param_table, sv_reg, pp_reg):
    """Calculate maximum C/N and C/P ratios for structural material.

    These ratios limit decomposition of structural material (i.e., material
    containing lignin). Lines 31-77 Predec.f

    Parameters:
        site_index_path (string): path to site spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        sv_reg (dict): map of key, path pairs giving paths to state
            variables for the current month
        pp_reg (dict): map of key, path pairs giving paths to persistent
            intermediate parameters that do not change over the course of
            the simulation.

    Modifies the persistent parameter rasters indexed by the following
    keys:
        pp_reg['rnewas_1_1_path']
        pp_reg['rnewas_1_2_path']
        pp_reg['rnewas_2_1_path']
        pp_reg['rnewas_2_2_path']
        pp_reg['rnewbs_1_1_path']
        pp_reg['rnewbs_1_2_path']
        pp_reg['rnewbs_2_1_path']
        pp_reg['rnewbs_2_2_path']

    Returns:
        None
    """
    strucc_1_nodata = pygeoprocessing.get_raster_info(
            sv_reg['strucc_1_path'])['nodata'][0]

    # temporary parameter rasters for structural ratios calculations
    temp_dir = tempfile.mkdtemp()
    param_val_dict = {}
    for iel in [1, 2]:
        for val in[
                'pcemic1_2', 'pcemic1_1', 'pcemic1_3', 'pcemic2_2',
                'pcemic2_1', 'pcemic2_3', 'rad1p_1', 'rad1p_2',
                'rad1p_3', 'varat1_1', 'varat22_1']:
            target_path = os.path.join(temp_dir, '{}_{}.tif'.format(val, iel))
            param_val_dict['{}_{}'.format(val, iel)] = target_path
            site_to_val = dict(
                [(site_code, float(table['{}_{}'.format(val, iel)])) for
                    (site_code, table) in site_param_table.iteritems()])
            pygeoprocessing.reclassify_raster(
                (site_index_path, 1), site_to_val, target_path,
                gdal.GDT_Float32, _IC_NODATA)

    def calc_rnewas_som1(
            pcemic1_2, pcemic1_1, pcemic1_3, struce_1, strucc_1):
        """Calculate C/<iel> ratio for decomposition into som1.

        This ratio is calculated separately for each nutrient (i.e., N, P).
        When material decomposes into the surface active organic pool, the
        C/<iel> ratio of decomposing material must be smaller than or equal to
        this ratio.

        Parameters:
            pcemic1_2 (numpy.ndarray): parameter, minimum C/<iel> ratio for
                surface active organic pool
            pcemic1_1 (numpy.ndarray): parameter, maximum C/<iel> ratio for
                surface active organic pool
            pcemic1_3 (numpy.ndarray): parameter, mimimum <iel> content of
                decomposing aboveground material, above which the C/<iel>
                ratio of the surface microbes equals pcemic1_2
            struce_1 (numpy.ndarray): state variable, <iel> in surface
                structural material
            strucc_1 (numpy.ndarray): state variable, C in surface
                structural material

        Returns:
            rnewas1, required ratio for decomposition of structural material
            into som1 for one nutrient
        """
        valid_mask = (
            (pcemic1_2 != _IC_NODATA) &
            (pcemic1_1 != _IC_NODATA) &
            (pcemic1_3 != _IC_NODATA) &
            (struce_1 != struce_1_nodata) &
            (strucc_1 != strucc_1_nodata))
        cemicb1 = numpy.empty(strucc_1.shape, dtype=numpy.float32)
        cemicb1[:] = _TARGET_NODATA
        cemicb1[valid_mask] = (
            (pcemic1_2[valid_mask] - pcemic1_1[valid_mask]) /
            pcemic1_3[valid_mask])

        rnewas1 = _aboveground_ratio(
            struce_1, strucc_1, pcemic1_1, pcemic1_2, pcemic1_3, cemicb1)
        return rnewas1

    def calc_rnewas_som2(
            pcemic2_2, pcemic2_1, pcemic2_3, struce_1, strucc_1, rad1p_1,
            rad1p_2, rad1p_3, pcemic1_2, rnewas1):
        """Calculate C/<iel> ratio for decomposition into som2.

        This ratio is calculated separately for each nutrient (i.e., N, P).
        When material decomposes into the surface slow organic pool, the
        C/<iel> ratio of decomposing material must be smaller than or equal to
        this ratio. A portion of the ratio of material entering som1, the
        surface active pool, is also added to som2 and calculated here.

        Parameters:
            pcemic2_2 (numpy.ndarray): parameter, minimum C/<iel> ratio for
                surface slow organic pool
            pcemic2_1 (numpy.ndarray): parameter, maximum C/<iel> ratio for
                surface slow organic pool
            pcemic2_3 (numpy.ndarray): parameter, mimimum <iel> content of
                decomposing aboveground material, above which the C/<iel>
                ratio of the surface slow organic pool equals pcemic1_2
            struce_1 (numpy.ndarray): state variable, <iel> in surface
                structural material
            strucc_1 (numpy.ndarray): state variable, C in surface
                structural material
            rad1p_1 (numpy.ndarray): parameter, intercept of regression used
                to calculate addition of <iel> from surface active pool
            rad1p_2 (numpy.ndarray): parameter, slope of regression used
                to calculate addition of <iel> from surface active pool
            rad1p_3 (numpy.ndarray): parameter, minimum allowable C/<iel>
                used to calculate addition term for C/<iel> ratio of som2
                formed from surface active pool
            pcemic1_2 (numpy.ndarray): parameter, minimum C/<iel> ratio for
                surface active organic pool
            rnewas1 (numpy.ndarray): derived, C/<iel> ratio for decomposition
                into som1

        Returns:
            rnewas2, required ratio for decomposition of structural material
            into som2 for one nutrient
        """
        valid_mask = (
            (pcemic2_2 != _IC_NODATA) &
            (pcemic2_1 != _IC_NODATA) &
            (pcemic2_3 != _IC_NODATA) &
            (struce_1 != struce_1_nodata) &
            (strucc_1 != strucc_1_nodata) &
            (rad1p_1 != _IC_NODATA) &
            (rad1p_2 != _IC_NODATA) &
            (rad1p_3 != _IC_NODATA) &
            (pcemic1_2 != _IC_NODATA) &
            (rnewas1 != _TARGET_NODATA))
        cemicb2 = numpy.empty(strucc_1.shape, dtype=numpy.float32)
        cemicb2[:] = _TARGET_NODATA
        cemicb2[valid_mask] = (
            (pcemic2_2[valid_mask] - pcemic2_1[valid_mask]) /
            pcemic2_3[valid_mask])

        rnewas2 = _aboveground_ratio(
            struce_1, strucc_1, pcemic2_1, pcemic2_2, pcemic2_3, cemicb2)

        radds1 = numpy.empty(strucc_1.shape, dtype=numpy.float32)
        radds1[:] = _TARGET_NODATA
        radds1[valid_mask] = (
            rad1p_1[valid_mask] + rad1p_2[valid_mask] *
            (rnewas1[valid_mask] - pcemic1_2[valid_mask]))
        rnewas2[valid_mask] = rnewas1[valid_mask] + radds1[valid_mask]
        rnewas2[valid_mask] = numpy.maximum(
            rnewas2[valid_mask], rad1p_3[valid_mask])
        return rnewas2

    for iel in [1, 2]:
        # calculate rnewas_iel_1 - aboveground material to SOM1
        struce_1_nodata = pygeoprocessing.get_raster_info(
            sv_reg['struce_1_{}_path'.format(iel)])['nodata'][0]
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
            (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_varat1_1,
            pp_reg['rnewbs_{}_1_path'.format(iel)],
            gdal.GDT_Float32, _TARGET_NODATA)
        # calculate rnewbs_iel_2 - belowground material to SOM2
        # rnewbs(iel,2) = varat22(1,iel)
        site_to_varat22_1 = dict([
            (site_code, float(table['varat22_1_{}'.format(iel)])) for
            (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_varat22_1,
            pp_reg['rnewbs_{}_2_path'.format(iel)],
            gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _yearly_tasks(
        site_index_path, site_param_table, aligned_inputs, month_index,
        year_reg):
    """Calculate quantities that remain static for 12 months.

    These quantities are annual precipitation and annual atmospheric N
    deposition. Century also calculates non-symbiotic soil N fixation once
    yearly, but here those were moved to monthly tasks.
    Century uses precipitation in the future 12 months (prcgrw) to predict
    root:shoot ratios, but here we instead use the sum of monthly
    precipitation in 12 months including the current one, if data for 12
    future months are not available.
    Lines 79-82 Eachyr.f

    Parameters:
        site_index_path (string): path to site spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including monthly precipitation
        month_index (int): current monthly step, relative to 0 so that
            month_index=0 at first monthly time step
        year_reg (dict): map of key, path pairs giving paths to the annual
            precipitation and N deposition rasters

    Modifies the rasters year_reg['annual_precip_path'] and
        year_reg['baseNdep_path']

    Returns:
        None

    Raises:
        ValueError if less than 12 monthly precipitation rasters can be found
    """
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

    precip_nodata = set([])
    for precip_raster in annual_precip_rasters:
        precip_nodata.update(
            set([pygeoprocessing.get_raster_info(precip_raster)['nodata'][0]]))
    if len(precip_nodata) > 1:
        raise ValueError("Precipitation rasters include >1 nodata value")
    precip_nodata = list(precip_nodata)[0]

    raster_sum(
        annual_precip_rasters, precip_nodata, year_reg['annual_precip_path'],
        _TARGET_NODATA)

    # calculate base N deposition
    # intermediate parameter rasters for this operation
    temp_dir = tempfile.mkdtemp()
    param_val_dict = {}
    for val in['epnfa_1', 'epnfa_2']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)

    def calc_base_N_dep(epnfa_1, epnfa_2, prcann):
        """Calculate base annual atmospheric N deposition.

        Parameters:
            epnfa_1 (numpy.ndarray): parameter, intercept of regression
                predicting atmospheric N deposition from precipitation
            epnfa_2 (numpy.ndarray): parameter, slope of regression predicting
                atmospheric N deposition from precipitation
            prcann (numpy.ndarray): derived, annual precipitation

        Returns:
            baseNdep, annual atmospheric N deposition
        """
        baseNdep = numpy.empty(prcann.shape, dtype=numpy.float32)
        baseNdep[:] = 0.
        valid_mask = (
            (epnfa_1 != _IC_NODATA) &
            (epnfa_2 != _IC_NODATA) &
            (prcann != _TARGET_NODATA))
        baseNdep[valid_mask] = (
            epnfa_1[valid_mask] +
            (epnfa_2[valid_mask] * numpy.minimum(prcann[valid_mask], 80.)))
        baseNdep[baseNdep < 0] = 0.
        return baseNdep

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['epnfa_1'], param_val_dict['epnfa_2'],
            year_reg['annual_precip_path']]],
        calc_base_N_dep, year_reg['baseNdep_path'], gdal.GDT_Float32,
        _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _shortwave_radiation(template_raster, month, shwave_path):
    """Calculate shortwave radiation outside the atmosphere.

    Shortwave radiation outside the atmosphere is calculated according to
    Penman (1948), "Natural evaporation from open water, bare soil and grass",
    Proc. Roy. Soc. London. The latitude of each pixel is required to
    calculate radiation and is calculated as an intermediate step from the
    input `template_raster`. shwave.f

    Parameters:
        template_raster (string): path to a raster in geographic coordinates
            that is aligned with model inputs
        month (int): current month of the year, such that month=0 indicates
            January
        shwave_path (string): path to shortwave radiation raster

    Modifies the raster indicated by `shwave_path`

    Returns:
        None
    """
    def shwave(month):
        def _shwave(latitude):
            """Calculate shortwave radiation outside the atmosphere.

            Parameters:
                latitude (float): latitude of current site in degrees
                month (int): current month of the year, such that month=1
                    indicates January

            Returns:
                shwave, short wave solar radiation outside the atmosphere
            """
            # Julian date in middle of each month of the year
            jday_list = [
                16, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350]
            jday = jday_list[month - 1]
            transcof = 0.8

            # Convert latitude from degrees to radians
            rlatitude = latitude * (numpy.pi / 180.0)

            # short wave solar radiation on a clear day
            declin = 0.401426 * numpy.sin(6.283185 * (jday - 77.0) / 365.0)

            temp = 1.0 - (-numpy.tan(rlatitude) * numpy.tan(declin))**2
            temp[temp < 0.] = 0.

            par1 = numpy.sqrt(temp)
            par2 = (-numpy.tan(rlatitude) * numpy.tan(declin))

            ahou = numpy.arctan2(par1, par2)
            ahou[ahou < 0.] = 0.

            solrad = (
                917.0 * transcof * (
                    ahou * numpy.sin(rlatitude) * numpy.sin(declin) +
                    numpy.cos(rlatitude) *
                    numpy.cos(declin) * numpy.sin(ahou)))

            # short wave radiation outside the atmosphere
            shwave = solrad / transcof
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

        target_band.WriteArray(
            y_vector, xoff=offset_map['xoff'], yoff=offset_map['yoff'])

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


def _reference_evapotranspiration(
        max_temp_path, min_temp_path, shwave_path, fwloss_4_path,
        pevap_path):
    """Calculate reference evapotranspiration.

    Reference evapotranspiration from the FAO Penman-Monteith equation in
    "Guidelines for computing crop water requirements", FAO Irrigation and
    drainage paper 56 (http://www.fao.org/docrep/X0490E/x0490e08.htm),
    modified by the parameter fwloss(4).

    Parameters:
            max_temp_path (string): path to maximum monthly temperature
            min_temp_path (string): path to minimum monthly temperature
            shwave_path (string): path to shortwave radiation outside the
                atmosphere
            fwloss_4_path (string): path to parameter, scaling factor for
                reference evapotranspiration
            pevap_path (string): path to result, reference evapotranspiration
                raster

    Modifies:
        The raster indicated by `pevap_path`

    Returns:
        None
    """
    def _calc_pevap(max_temp, min_temp, shwave, fwloss_4):
        """Calculate reference evapotranspiration.

        Pevap.f

        Parameters:
            max_temp (numpy.ndarray): input, maximum monthly temperature
            min_temp (numpy.ndarray): input, minimum monthly temperature
            shwave (numpy.ndarray): derived, shortwave radiation outside the
                atmosphere
            fwloss_4 (numpy.ndarray): parameter, scaling factor for reference
                evapotranspiration

        Returns:
            pevap, reference evapotranspiration
        """
        const1 = 0.0023
        const2 = 17.8
        langleys2watts = 54.0

        valid_mask = (
            (max_temp != maxtmp_nodata) &
            (min_temp != mintmp_nodata) &
            (shwave != _TARGET_NODATA) &
            (fwloss_4 != _IC_NODATA))
        trange = numpy.empty(fwloss_4.shape, dtype=numpy.float32)
        trange[:] = _TARGET_NODATA
        trange[valid_mask] = max_temp[valid_mask] - min_temp[valid_mask]
        tmean = numpy.empty(fwloss_4.shape, dtype=numpy.float32)
        tmean[:] = _IC_NODATA
        tmean[valid_mask] = (max_temp[valid_mask] + min_temp[valid_mask]) / 2.0

        # daily reference evapotranspiration
        daypet = numpy.empty(fwloss_4.shape, dtype=numpy.float32)
        daypet[:] = _TARGET_NODATA
        in1 = const1 * (tmean[valid_mask] + const2)
        in2 = numpy.sqrt(trange[valid_mask])
        in3 = (shwave[valid_mask] / langleys2watts)
        daypet[valid_mask] = (
            const1 * (tmean[valid_mask] + const2) *
            numpy.sqrt(trange[valid_mask]) *
            (shwave[valid_mask] / langleys2watts))

        # monthly reference evapotranspiration, from mm to cm,
        # bounded to be at least 0.5
        monpet = (daypet * 30.) / 10.
        monpet[monpet <= 0.5] = 0.5

        pevap = numpy.empty(fwloss_4.shape, dtype=numpy.float32)
        pevap[:] = _TARGET_NODATA
        pevap[valid_mask] = monpet[valid_mask] * fwloss_4[valid_mask]
        return pevap

    maxtmp_nodata = pygeoprocessing.get_raster_info(
        max_temp_path)['nodata'][0]
    mintmp_nodata = pygeoprocessing.get_raster_info(
        min_temp_path)['nodata'][0]
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            max_temp_path, min_temp_path, shwave_path, fwloss_4_path]],
        _calc_pevap, pevap_path, gdal.GDT_Float32, _TARGET_NODATA)


def _potential_production(
        aligned_inputs, site_param_table, current_month, month_index,
        pft_id_set, veg_trait_table, sv_reg, pp_reg, month_reg):
    """Calculate above- and belowground potential production.

    Potential production of each plant functional type is calculated
    as total potential production given incoming solar radiation,
    limited by temperature, soil moisture, and obstruction by biomass and
    litter. Further modification of potential production according to
    limitation by water and nutrient availability is calculated in the
    root:shoot ratio submodel. Lines 57-148 Potcrp.f

    Parameters:
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including precipitation, temperature,
            plant functional type composition, and site spatial index
        site_param_table (dict): map of site spatial indices to dictionaries
            containing site parameters
        current_month (int): month of the year, such that current_month=1
            indicates January
        month_index (int): month of the simulation, such that month_index=13
            indicates month 13 of the simulation
        pft_id_set (set): set of integers identifying plant functional types
        veg_trait_table (dict): map of pft id to dictionaries containing
            plant functional type parameters
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month
        pp_reg (dict): map of key, path pairs giving paths to persistent
            intermediate parameters that do not change over the course of
            the simulation
        month_reg (dict): map of key, path pairs giving paths to intermediate
            calculated values that are shared between submodels

    Modifies:
        The raster indicated by `month_reg['h2ogef_1_<PFT>']` for each
            plant functional type (PFT) where growth is scheduled to occur in
            this month
        The raster indicated by `month_reg['tgprod_pot_prod_<PFT>']` for each
            plant functional type (PFT) where growth is scheduled to occur in
            this month

    Returns:
        None
    """
    # if growth does not occur this month for all PFTs,
    # skip the rest of the function
    do_PFT = []
    for pft_index in pft_id_set:
        if str(current_month) in veg_trait_table[pft_index]['growth_months']:
            do_PFT.append(pft_index)
    if not do_PFT:
        return

    def multiply_positive_rasters(raster1, raster2):
        """Multiply positive values in two rasters."""
        valid_mask = (raster1 > 0.) & (raster2 > 0.)
        result = numpy.zeros(raster1.shape)
        result[valid_mask] = raster1[valid_mask] * raster2[valid_mask]
        return result

    def calc_ctemp(aglivc, pmxbio, maxtmp, pmxtmp, mintmp, pmntmp):
        """Calculate soil temperature relative to its effect on growth.

        Soil temperature is calculated from monthly temperature inputs and
        modified by total standing live biomass. Lines 69-84 Potcrp.f

        Parameters:
            aglivc (numpy.ndarray): derived, sum of aglivc (carbon in
                aboveground live biomass) across plant functional types
            pmxbio (numpy.ndarray): parameter, maximum biomass impact on
                temperature
            maxtmp (numpy.ndarray): derived, average maximum monthly
                temperature
            pmxtmp (numpy.ndarray): parameter, scaling factor for effect of
                biomass on monthly maximum temperature
            mintmp (numpy.ndarray): input, average minimum monthly temperature
            pmntmp (numpy.ndarray): parameter, scaling factor for effect of
                biomass on monthly minimum temperature

        Returns:
            ctemp, effect of soil temperature on potential production

        """
        bio = numpy.empty(aglivc.shape, dtype=numpy.float32)
        bio[:] = _IC_NODATA
        valid_mask = (
            (aglivc >= 0.) &
            (pmxbio != _IC_NODATA) &
            (maxtmp != maxtmp_nodata) &
            (pmxtmp != _IC_NODATA) &
            (mintmp != mintmp_nodata) &
            (pmntmp != _IC_NODATA))
        bio[valid_mask] = aglivc[valid_mask] * 2.5
        bio[bio > pmxbio] = pmxbio[bio > pmxbio]
        bio[pmxbio < 0] = _IC_NODATA

        # Maximum temperature
        tmxs = numpy.empty(aglivc.shape, dtype=numpy.float32)
        tmxs[:] = _IC_NODATA
        tmxs[valid_mask] = (
            maxtmp[valid_mask] + (
                (25.4/(1. + 18. * numpy.exp(-0.20 * maxtmp[valid_mask]))) *
                (numpy.exp(pmxtmp[valid_mask] * bio[valid_mask]) - 0.13)))

        # Minimum temperature
        tmns = numpy.empty(aglivc.shape, dtype=numpy.float32)
        tmns[:] = _IC_NODATA
        tmns[valid_mask] = (
            mintmp[valid_mask] +
            (pmntmp[valid_mask] * bio[valid_mask] - 1.78))

        # Average temperature
        ctemp = numpy.empty(aglivc.shape, dtype=numpy.float32)
        ctemp[:] = _IC_NODATA
        ctemp[valid_mask] = (tmxs[valid_mask] + tmns[valid_mask])/2.
        return ctemp

    def calc_potprd(ctemp, ppdf_1, ppdf_2, ppdf_3, ppdf_4):
        """Calculate the limiting effect of temperature on growth.

        Estimated soil temperature restricts potential production according to
        a Poisson Density Function curve described by the plant functional
        type-specific parameters ppdf_1-4.. Lines 73-84 Potcrp.f

        Parameters:
            ctemp (numpy.ndarray): derived, soil temperature as calculated from
                monthly temperature and modified by standing live biomass
            ppdf_1 (numpy.ndarray): parameter, optimum temperature for growth
            ppdf_2 (numpy.ndarray): parameter, maximum temperature for growth
            ppdf_3 (numpy.ndarray): parameter, left curve shape for Poisson
                Density Function curve describing growth as function of
                temperature
            ppdf_4 (numpy.ndarray): parameter, right curve shape for Poisson
                Density Function curve describing growth as function of
                temperature

        Returns:
            potprd, scaling factor describing potential production limited
                by temperature
        """
        valid_mask = (
            (ctemp != _IC_NODATA) &
            (ppdf_1 != _IC_NODATA) &
            (ppdf_2 != _IC_NODATA) &
            (ppdf_3 != _IC_NODATA) &
            (ppdf_4 != _IC_NODATA))
        frac = numpy.empty(ctemp.shape, dtype=numpy.float32)
        frac[:] = _TARGET_NODATA
        frac[valid_mask] = (
            (ppdf_2[valid_mask] - ctemp[valid_mask]) /
            (ppdf_2[valid_mask] - ppdf_1[valid_mask]))
        gpdf = numpy.empty(ctemp.shape, dtype=numpy.float32)
        gpdf[:] = _TARGET_NODATA
        gpdf[valid_mask] = (numpy.exp(
            (ppdf_3[valid_mask]/ppdf_4[valid_mask]) *
            (1. - numpy.power(frac[valid_mask], ppdf_4[valid_mask]))) *
            numpy.power(frac[valid_mask], ppdf_3[valid_mask]))
        return gpdf

    def calc_h2ogef_1(
            pevap, avh2o_1, precip, wc, pprpts_1, pprpts_2, pprpts_3):
        """Calculate the limiting factor of water availability on growth.

        Soil moisture restricts potential production according to the ratio
        of available water to reference evapotranspiration. The shape of the
        linear relationship of this ratio to potential production is
        controlled by the site parameters pprpts_1, pprpts_2, and pprpts_3.
        Lines 57-64 Potcrp.f

        Parameters:
            pevap (numpy.ndarray): derived, reference evapotranspiration
            avh2o_1 (numpy.ndarray): derived, water available to this plant
                functional type for growth
            precip (numpy.ndarray): input, precipitation for the current month
            wc (numpy.ndarray): derived, water content in soil layer 1
            pprpts_1 (numpy.ndarray): parameter, the minimum ratio of
                available water to reference evapotranspiration that limits
                production completely
            pprpts_2 (numpy.ndarray): parameter, influences the slope of the
                line predicting potential production from available water
            pprpts_3 (numpy.ndarray): parameter, the ratio of available water
                to reference evapotranspiration above which production is
                not restricted

        Returns:
            h2ogef_1, scaling factor describing potential production limited
                by soil moisture

        """
        valid_mask = (
            (pevap != _TARGET_NODATA) &
            (avh2o_1 != avh2o_1_nodata) &
            (precip != precip_nodata) &
            (wc != _TARGET_NODATA) &
            (pprpts_1 != _IC_NODATA) &
            (pprpts_2 != _IC_NODATA) &
            (pprpts_3 != _IC_NODATA))
        h2ogef_prior = numpy.empty(pevap.shape, dtype=numpy.float32)
        h2ogef_prior[:] = _TARGET_NODATA
        h2ogef_prior[valid_mask] = numpy.where(
            pevap[valid_mask] >= 0.01,
            (avh2o_1[valid_mask] + precip[valid_mask])/pevap[valid_mask],
            0.01)

        intcpt = pprpts_1 + (pprpts_2 * wc)
        slope = 1. / (pprpts_3 - intcpt)

        h2ogef_1 = numpy.empty(pevap.shape, dtype=numpy.float32)
        h2ogef_1[:] = _TARGET_NODATA
        h2ogef_1[valid_mask] = (
            1.0 + slope[valid_mask] *
            (h2ogef_prior[valid_mask] - pprpts_3[valid_mask]))

        h2ogef_1[valid_mask] = numpy.clip(h2ogef_1[valid_mask], 0.01, 1.)
        return h2ogef_1

    def calc_biof(sum_stdedc, sum_aglivc, strucc_1, pmxbio, biok5):
        """Calculate the effect of obstruction on growth.

        Live biomass, standing dead biomass, and litter reduce potential
        production through obstruction. The shape of the relationship between
        standing biomass and litter and potential production is controlled by
        the site parameter pmxbio and the plant functional type parameter
        biok5. Lines 91-120 Potcrp.f

        Parameters:
            sum_stdedc (numpy.ndarray): derived, total carbon in standing dead
                biomass across plant functional types
            sum_aglivc (numpy.ndarray): derived, total carbon in aboveground
                live biomass across plant functional types
            strucc_1 (numpy.ndarray): derived, carbon in surface litter
            pmxbio (numpy.ndarray): parameter, maximum biomass impact on
                potential production
            biok5 (numpy.ndarray): parameter, level of standing dead biomass
                and litter

        Returns:
            biof, scaling factor describing potential production limited
                by obstruction
        """
        valid_mask = (
            (strucc_1 != strucc_1_nodata) &
            (pmxbio != _IC_NODATA) &
            (biok5 != _IC_NODATA))

        bioc = numpy.empty(sum_stdedc.shape, dtype=numpy.float32)
        bioc[:] = _IC_NODATA
        bioc[valid_mask] = numpy.where(
            ((sum_stdedc[valid_mask] + 0.1*strucc_1[valid_mask]) <= 0.), 0.01,
            (sum_stdedc[valid_mask] + 0.1*strucc_1[valid_mask]))
        bioc[valid_mask] = numpy.where(
            (bioc[valid_mask] > pmxbio[valid_mask]), pmxbio[valid_mask],
            bioc[valid_mask])

        bioprd = numpy.empty(sum_stdedc.shape, dtype=numpy.float32)
        bioprd[:] = _IC_NODATA
        bioprd[valid_mask] = 1. - (
            bioc[valid_mask] / (biok5[valid_mask] + bioc[valid_mask]))

        temp1 = 1. - bioprd
        temp2 = temp1 * 0.75
        temp3 = temp1 * 0.25

        ratlc = numpy.empty(sum_stdedc.shape, dtype=numpy.float32)
        ratlc[:] = _IC_NODATA
        ratlc[valid_mask] = sum_aglivc[valid_mask] / bioc[valid_mask]

        biof = numpy.empty(sum_stdedc.shape, dtype=numpy.float32)
        biof[:] = _TARGET_NODATA
        biof[valid_mask] = numpy.where(
            ratlc[valid_mask] <= 1.,
            (bioprd[valid_mask] + (temp2[valid_mask] * ratlc[valid_mask])),
            numpy.where(
                ratlc[valid_mask] <= 2.,
                (bioprd[valid_mask] + temp2[valid_mask]) +
                temp3[valid_mask] * (ratlc[valid_mask] - 1.),
                1.))
        return biof

    def calc_tgprod_pot_prod(prdx_1, shwave, potprd, h2ogef_1, biof):
        """Calculate total potential production.

        Total above- and belowground potential biomass production is calculated
        as the total potential production given solar radiation and the
        intrinsinc growth capacity of the plant functional type, modified by
        limiting factors of temperature, soil moisture, and obstruction by
        standing biomass and litter. Line 147 Potcrp.f

        Parameters:
            prdx_1 (numpy.ndarray): parameter, the intrinsic capacity of the
                plant functional type for growth per unit of solar radiation
            shwave (numpy.ndarray): derived, shortwave solar radiation outside
                the atmosphere
            potprd (numpy.ndarray): parameter, scaling factor describing
                limiting effect of temperature
            h2ogef_1 (numpy.ndarray): derived, scaling factor describing the
                limiting effect of soil moisture
            biof (numpy.ndarray): derived, scaling factor describing the
                limiting effect of obstruction by standing biomass and litter

        Returns:
            tgprod_pot_prod, total above- and belowground potential biomass
                production
        """
        valid_mask = (
            (prdx_1 != _IC_NODATA) &
            (shwave != _TARGET_NODATA) &
            (potprd != _TARGET_NODATA) &
            (h2ogef_1 != _TARGET_NODATA) &
            (biof != _TARGET_NODATA))
        tgprod_pot_prod = numpy.empty(prdx_1.shape, dtype=numpy.float32)
        tgprod_pot_prod[:] = _TARGET_NODATA

        tgprod_pot_prod[valid_mask] = (
            prdx_1[valid_mask] * shwave[valid_mask] * potprd[valid_mask] *
            h2ogef_1[valid_mask] * biof[valid_mask])
        return tgprod_pot_prod

    # temporary intermediate rasters for calculating total potential production
    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {}
    # site-level temporary calculated values
    for val in ['sum_aglivc', 'sum_stdedc', 'ctemp', 'shwave', 'pevap']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))
    # PFT-level temporary calculated values
    for pft_i in pft_id_set:
        for val in [
                'aglivc_weighted', 'stdedc_weighted', 'potprd', 'biof']:
            temp_val_dict['{}_{}'.format(val, pft_i)] = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))

    # temporary parameter rasters for calculating total potential production
    param_val_dict = {}
    # site-level parameters
    for val in [
            'pmxbio', 'pmxtmp', 'pmntmp', 'fwloss_4', 'pprpts_1',
            'pprpts_2', 'pprpts_3']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (aligned_inputs['site_index'], 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)
    # PFT-level parameters
    for val in [
            'ppdf_1', 'ppdf_2', 'ppdf_3', 'ppdf_4', 'biok5', 'prdx_1']:
        for pft_i in do_PFT:
            target_path = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
            param_val_dict['{}_{}'.format(val, pft_i)] = target_path
            fill_val = veg_trait_table[pft_i][val]
            pygeoprocessing.new_raster_from_base(
                aligned_inputs['site_index'], target_path, gdal.GDT_Float32,
                [_IC_NODATA], fill_value_list=[fill_val])

    # calculate intermediate quantities that do not differ between PFTs:
    # sum of aglivc (standing live biomass) and stdedc (standing dead biomass)
    # across PFTs, weighted by % cover of each PFT
    maxtmp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['max_temp_{}'.format(current_month)])['nodata'][0]
    mintmp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['min_temp_{}'.format(current_month)])['nodata'][0]
    precip_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['precip_{}'.format(month_index)])['nodata'][0]
    strucc_1_nodata = pygeoprocessing.get_raster_info(
        sv_reg['strucc_1_path'])['nodata'][0]
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

        raster_sum(
            weighted_path_list, _TARGET_NODATA,
            temp_val_dict['sum_{}'.format(sv)], _TARGET_NODATA,
            nodata_remove=True)

    # ctemp, soil temperature relative to impacts on growth
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['sum_aglivc'],
            param_val_dict['pmxbio'],
            aligned_inputs['max_temp_{}'.format(current_month)],
            param_val_dict['pmxtmp'],
            aligned_inputs['min_temp_{}'.format(current_month)],
            param_val_dict['pmntmp']]],
        calc_ctemp, temp_val_dict['ctemp'], gdal.GDT_Float32, _IC_NODATA)

    # shwave, shortwave radiation outside the atmosphere
    _shortwave_radiation(
        aligned_inputs['site_index'], current_month, temp_val_dict['shwave'])

    # pet, reference evapotranspiration modified by fwloss parameter
    _reference_evapotranspiration(
        aligned_inputs['max_temp_{}'.format(current_month)],
        aligned_inputs['min_temp_{}'.format(current_month)],
        temp_val_dict['shwave'],
        param_val_dict['fwloss_4'],
        temp_val_dict['pevap'])

    # calculate quantities that differ between PFTs
    for pft_i in do_PFT:
        avh2o_1_nodata = pygeoprocessing.get_raster_info(
            sv_reg['avh2o_1_{}_path'.format(pft_i)])['nodata'][0]
        # potprd, the limiting effect of temperature
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in
                temp_val_dict['ctemp'],
                param_val_dict['ppdf_1_{}'.format(pft_i)],
                param_val_dict['ppdf_2_{}'.format(pft_i)],
                param_val_dict['ppdf_3_{}'.format(pft_i)],
                param_val_dict['ppdf_4_{}'.format(pft_i)]],
            calc_potprd, temp_val_dict['potprd_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

        # h2ogef_1, the limiting effect of soil water availability
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in
                temp_val_dict['pevap'],
                sv_reg['avh2o_1_{}_path'.format(pft_i)],
                aligned_inputs['precip_{}'.format(month_index)],
                pp_reg['wc_path'],
                param_val_dict['pprpts_1'],
                param_val_dict['pprpts_2'],
                param_val_dict['pprpts_3']],
            calc_h2ogef_1, month_reg['h2ogef_1_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

        # biof, the limiting effect of obstruction
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in
                temp_val_dict['sum_stdedc'],
                temp_val_dict['sum_aglivc'],
                sv_reg['strucc_1_path'],
                param_val_dict['pmxbio'],
                param_val_dict['biok5_{}'.format(pft_i)]],
            calc_biof, temp_val_dict['biof_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

        # total potential production
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in
                param_val_dict['prdx_1_{}'.format(pft_i)],
                temp_val_dict['shwave'],
                temp_val_dict['potprd_{}'.format(pft_i)],
                month_reg['h2ogef_1_{}'.format(pft_i)],
                temp_val_dict['biof_{}'.format(pft_i)]],
            calc_tgprod_pot_prod,
            month_reg['tgprod_pot_prod_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _calc_favail_P(sv_reg, param_val_dict):
    """Calculate the fraction of P in surface layer available to plants.

    This must be performed after the sum of mineral N in the surface layer
    is calculated because the fraction of labile P available to plants is
    impacted by the amount of mineral N in the surface layer.

    Parameters:
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month, including minerl_1_1, mineral N in the
            surface layer
        param_val_dict (dict): map of key, path pairs giving paths to
            site-level parameters, including favail_4, favail_5, favail_6,
            and favail_2

    Modifies:
        The raster indicated by `param_val_dict['favail_2']`

    Returns:
        None
    """
    def favail_P_op(minerl_1_1, favail_4, favail_5, favail_6):
        """Calculate the fraction of P in surface layer available to plants.

        The fraction of labile P available to plants depends on mineral N in
        the surface layer and site parameters favail_4, favail_5, favail_6.
        Line 395 Simsom.f

        Parameters:
            minerl_1_1 (numpy.ndarray): state variable, mineral N in the
                surface layer
            favail_4 (numpy.ndarray): parameter, minimum fraction of P
                available
            favail_5 (numpy.ndarray): parameter, maximum fraction of P
                available
            favail_6 (numpy.ndarray): parameter, mineral N in surface layer
                required to attain maximum fraction of P available

        Returns:
            favail_P, fraction of mineral P available to plants
        """
        valid_mask = (
            (minerl_1_1 != minerl_1_1_nodata) &
            (favail_4 != _IC_NODATA) &
            (favail_5 != _IC_NODATA) &
            (favail_6 != _IC_NODATA))

        interim = numpy.empty(minerl_1_1.shape, dtype=numpy.float32)
        interim[:] = _IC_NODATA

        interim[valid_mask] = (
            favail_4[valid_mask] + minerl_1_1[valid_mask] *
            (favail_5[valid_mask] - favail_4[valid_mask]) /
            favail_6[valid_mask])

        favail_P = numpy.empty(minerl_1_1.shape, dtype=numpy.float32)
        favail_P[:] = _IC_NODATA

        favail_P[valid_mask] = numpy.maximum(
            favail_4[valid_mask], numpy.minimum(
                interim[valid_mask], favail_5[valid_mask]))
        return favail_P

    minerl_1_1_nodata = pygeoprocessing.get_raster_info(
        sv_reg['minerl_1_1_path'])['nodata'][0]
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in
            sv_reg['minerl_1_1_path'],
            param_val_dict['favail_4'],
            param_val_dict['favail_5'],
            param_val_dict['favail_6']],
        favail_P_op, param_val_dict['favail_2'],
        gdal.GDT_Float32, _IC_NODATA)


def _calc_available_nutrient(
        pft_i, iel, pft_param_dict, sv_reg, site_param_table, site_index_path,
        favail_path, tgprod_path, eavail_path):
    """Calculate nutrient available to a plant functional type.

    The nutrient available is the sum of mineral nutrient (N or P) in soil
    layers accessible by the roots of the plant functional type, modified
    by the fraction of nutrient available to plants and the current root
    biomass.

    Parameters:
        pft_i (int): plant functional type index
        iel (int): nutrient index (iel=1 indicates N, iel=2 indicates P)
        pft_param_dict (dict): map of key, value pairs giving the values of
            parameters for this plant functional type
            (i.e., veg_trait_table[pft_i] for this pft_i)
        sv_reg (dict): map of key, path pairs giving paths to state
            variables for the current month
        site_index_path (string): path to site spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        favail_path (string): path to raster containing the appropriate value
            of the parameter favail. For nitrogen, this parameter is supplied
            directly as user input, but for phosphorus, it must be calculated
            from other parameters.
        tgprod_path (string): path to raster containing total potential
            production limited by shortwave radiation, temperature, and
            moisture
        eavail_path (string): path to location to store the result, nutrient
            available to the plant functional type

    Modifies:
        the raster indicated by `eavail_path`

    Returns:
        None
    """
    def calc_eavail(rictrl, bglivc, riint, availm, favail, crpstg):
        """Calculate available nutrient.

        Parameters:
            rictrl (numpy.ndarray): parameter, scaling factor used to
                calculate the impact of root biomass on nutrient availability
            bglivc (numpy.ndarray): state variable, carbon in belowground
                live biomass
            riint (numpy.ndarray): parameter, intercept used to calculate the
                impact of root biomass on nutrient availability
            availm (numpy.ndarray): derived, the sum of mineral nutrient in
                soil layers accessible by this plant functional type
            favail (numpy.ndarray): parameter, fraction of the nutrient
                available each month to plants
            crpstg (numpy.ndarray): state variable, nutrient in
                retranslocation storage pool for the plant functional type

        Returns:
            eavail, the nutrient available to the plant functional type
        """
        valid_mask = (
            (rictrl != _IC_NODATA) &
            (bglivc != bglivc_nodata) &
            (riint != _IC_NODATA) &
            (favail != _IC_NODATA) &
            (crpstg != crpstg_nodata))

        rimpct = numpy.empty(rictrl.shape, dtype=numpy.float32)
        rimpct[:] = _TARGET_NODATA
        rimpct[valid_mask] = numpy.where(
            ((rictrl[valid_mask] * bglivc[valid_mask] * 2.5) > 33.),
            1., 1. - riint[valid_mask] * numpy.exp(
                -rictrl[valid_mask] * bglivc[valid_mask] * 2.5))

        eavail = numpy.empty(rictrl.shape, dtype=numpy.float32)
        eavail[:] = _TARGET_NODATA
        eavail[valid_mask] = (
            (availm[valid_mask] * favail[valid_mask] * rimpct[valid_mask]) +
            crpstg[valid_mask])
        return eavail

    def add_symbiotic_fixed_N(eavail_prior, snfxmx, tgprod):
        """Add nitrogen fixed by the plant to nutrient available.

        Some nitrogen may be fixed by the plant, and this must be added
        to available mineral nitrogen. Nitrogen fixed by the plant is
        calculated from total potential production and the maximum
        rate of N fixation.

        Parameters:
            eavail_prior (numpy.ndarray): derived, mineral nitrogen available
                to the plant functional type, calculated with calc_eavail()
            snfxmx (numpy.ndarray): parameter, maximum rate of symbiotic
                nitrogen fixation
            tgprod (numpy.ndarray): derived, total above- and belowground
                potential production

        Returns:
            eavail, total N available including N fixed by the plant
        """
        valid_mask = (
            (eavail_prior != _TARGET_NODATA)  &
            (snfxmx != _IC_NODATA)  &
            (tgprod != _TARGET_NODATA))

        maxNfix = numpy.empty(eavail_prior.shape, dtype=numpy.float32)
        maxNfix[:] = _TARGET_NODATA
        maxNfix[valid_mask] = snfxmx[valid_mask] * (tgprod[valid_mask] / 2.5)

        eavail = numpy.empty(eavail_prior.shape, dtype=numpy.float32)
        eavail[:] = _TARGET_NODATA
        eavail[valid_mask] = eavail_prior[valid_mask] + maxNfix[valid_mask]
        return eavail

    bglivc_nodata = pygeoprocessing.get_raster_info(
        sv_reg['bglivc_{}_path'.format(pft_i)])['nodata'][0]
    crpstg_nodata = pygeoprocessing.get_raster_info(
        sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)])['nodata'][0]

    # temporary intermediate rasters for calculating available nutrient
    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {}
    for val in ['availm']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))

    param_val_dict = {}
    for val in ['rictrl', 'riint']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)
    for val in ['snfxmx_1']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        fill_val = pft_param_dict[val]
        pygeoprocessing.new_raster_from_base(
            site_index_path, target_path, gdal.GDT_Float32,
            [_IC_NODATA], fill_value_list=[fill_val])

    nlay = int(pft_param_dict['nlaypg'])
    mineral_raster_list = [
        sv_reg['minerl_{}_{}_path'.format(lyr, iel)] for lyr in xrange(
            1, nlay + 1)]

    mineral_nodata = set([])
    for mineral_raster in mineral_raster_list:
        mineral_nodata.update(
            set([pygeoprocessing.get_raster_info(
                mineral_raster)['nodata'][0]]))
    if len(mineral_nodata) > 1:
        raise ValueError(
            "Mineral rasters for element {} contain >1 nodata value".format(
                iel))
    raster_sum(
        mineral_raster_list, mineral_nodata, temp_val_dict['availm'],
        _TARGET_NODATA, nodata_remove=True)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['rictrl'],
            sv_reg['bglivc_{}_path'.format(pft_i)],
            param_val_dict['riint'],
            temp_val_dict['availm'],
            favail_path,
            sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)]]],
        calc_eavail, eavail_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    if iel == 1:
        eavail_prior_path = os.path.join(temp_dir, 'eavail_prior.tif')
        shutil.copyfile(eavail_path, eavail_prior_path)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                eavail_prior_path,
                param_val_dict['snfxmx_1'],
                tgprod_path]],
            add_symbiotic_fixed_N, eavail_path,
            gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _calc_nutrient_demand(
        biomass_production_path, fraction_allocated_to_roots_path,
        cercrp_min_above_path, cercrp_min_below_path, demand_path):
    """Calculate the demand of one nutrient by a plant functional type.

    Demand is calculated from total biomass production, the fraction of biomass
    production allocated to roots, and the minimum carbon/nutrient ratios of
    above- and belowground live biomass.  Lines 88-92 CropDynC.f and line
    65, Nutrlm.f

    Parameters:
        biomass_production_path (string): path to raster giving total
            biomass production
        fraction_allocated_to_roots_path (string): path to raster giving
            the fraction fo total biomass production allocated to roots
        cercrp_min_above_path (string): path to raster giving the minimum
            ratio of carbon to nutrient in aboveground live biomass
        cercrp_min_below_path (string): path to raster giving the minimum
            ratio of carbon to nutrient in belowground live biomass

    Modifies:
        The raster indicated by `demand_path`

    Returns:
        None
    """
    def nutrient_demand_op(
            biomass_production, root_fraction, cercrp_min_above,
            cercrp_min_below):
        """Calculate nutrient demand.

        Parameters:
            biomass_production (numpy.ndarray): derived, total biomass
                production
            root_fraction (numpy.ndarray): derived, fraction of biomass
                allocated to roots
            cercrp_min_above (numpy.ndarray): derived, minimum carbon to
                nutrient ratio of new aboveground live material
            cercrp_min_below (numpy.ndarray): derived, minimum carbon to
                nutrient ratio of new belowground live material

        Returns:
            demand_e, nutrient demand
        """
        valid_mask = (
            (biomass_production != _TARGET_NODATA) &
            (root_fraction != _TARGET_NODATA) &
            (cercrp_min_above != _TARGET_NODATA) &
            (cercrp_min_below != _TARGET_NODATA))

        demand_above = numpy.empty(root_fraction.shape, dtype=numpy.float32)
        demand_above[:] = _TARGET_NODATA
        demand_above[valid_mask] = (
            ((biomass_production[valid_mask] *
                (1. - root_fraction[valid_mask])) / 2.5) *
            (1. / cercrp_min_above[valid_mask]))

        demand_below = numpy.empty(root_fraction.shape, dtype=numpy.float32)
        demand_below[:] = _TARGET_NODATA
        demand_below[valid_mask] = (
            ((biomass_production[valid_mask] *
                (root_fraction[valid_mask])) / 2.5) *
            (1. / cercrp_min_below[valid_mask]))

        demand_e = numpy.empty(root_fraction.shape, dtype=numpy.float32)
        demand_e[:] = _TARGET_NODATA
        demand_e[valid_mask] = (
            demand_above[valid_mask] + demand_below[valid_mask])
        return demand_e

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in
            biomass_production_path, fraction_allocated_to_roots_path,
            cercrp_min_above_path, cercrp_min_below_path],
        nutrient_demand_op, demand_path,
        gdal.GDT_Float32, _TARGET_NODATA)


def calc_provisional_fracrc(
        annual_precip, frtcindx, bgppa, bgppb, agppa, agppb,
        cfrtcw_1, cfrtcw_2, cfrtcn_1, cfrtcn_2):
    """Calculate provisional fraction of carbon allocated to roots.

    A temporary provisional fraction of carbon allocated to roots must be
    calculated prior to calculating plant demand for N and P. The value
    of this provisional fraction depends on whether the plant functional
    type is modeled as a perennial plant or with the "Great Plains"
    equation of Parton et al. 1987, "Analysis of factors controlling soil
    organic matter levels in Great Plains grasslands", Soil Science
    Society of America Journal. Lines 36-47 cropDynC.f

    Parameters:
        annual_precip (numpy.ndarray): derived, sum of monthly
            precipitation over twelve months including the current month
        frtcindx (numpy.ndarray): parameter, flag indicating whether
            root:shoot allocation follows the Great Plains equation
            (frtcindx=0) or as a perennial plant (frtcindx=1)
        bgppa (numpy.ndarray): parameter, intercept in regression
            estimating belowground production from annual precipitation
            if frtcindx=0
        bgppb (numpy.ndarray): parameter, slope in regression estimating
            belowground production from annual precipitation if
            frtcindx=0
        agppa (numpy.ndarray): parameter, intercept in regression
            estimating aboveground production from annual precipitation
            if frtcindx=0
        agppb (numpy.ndarray): parameter, slope in regression estimating
            aboveground production from annual precipitation if
            frtcindx=0
        cfrtcw_1 (numpy.ndarray): parameter, maximum fraction of carbon
            allocated to roots under maximum water stress if frtcindx=1
        cfrtcw_2 (numpy.ndarray): parameter, minimum fraction of carbon
            allocated to roots without water stress if frtcindx=1
        cfrtcn_1 (numpy.ndarray): parameter, maximum fraction of carbon
            allocated to roots under maximum nutrient stress if frtcindx=1
        cfrtcn_2 (numpy.ndarray): parameter, minimum fraction of carbon
            allocated to roots under no nutrient stress if frtcindx=1

    Returns:
        fracrc_p, provisional fraction of carbon allocated to roots
    """
    valid_mask = (
        (annual_precip != _TARGET_NODATA) &
        (frtcindx != _IC_NODATA))
    rtsh = numpy.empty(annual_precip.shape, dtype=numpy.float32)
    rtsh[:] = _TARGET_NODATA
    rtsh[valid_mask] = (
        (bgppa[valid_mask] +
            annual_precip[valid_mask] * bgppb[valid_mask]) /
        (agppa[valid_mask] + annual_precip[valid_mask] *
            agppb[valid_mask]))

    fracrc_p = numpy.empty(annual_precip.shape, dtype=numpy.float32)
    fracrc_p[:] = _TARGET_NODATA
    fracrc_p[valid_mask] = numpy.where(
        frtcindx[valid_mask] == 0,
        (1.0 / (1.0 / rtsh[valid_mask] + 1.0)),
        ((cfrtcw_1[valid_mask] + cfrtcw_2[valid_mask] +
            cfrtcn_1[valid_mask] + cfrtcn_2[valid_mask]) / 4.0))
    return fracrc_p


def calc_ce_ratios(
        pramn_1_path, pramn_2_path, aglivc_path, biomax_path,
        pramx_1_path, pramx_2_path, prbmn_1_path, prbmn_2_path,
        prbmx_1_path, prbmx_2_path, annual_precip_path, month_reg,
        pft_i, iel):
    """Calculate minimum and maximum carbon to nutrient ratios.

    Minimum and maximum C/E ratios are used to calculate demand for a
    nutrient by a plant functional type. This function calculates the
    ratios for above- and belowground plant portions, for one plant
    functional type and one nutrient. Fltce.f

    Parameters:
        pramn_1_path (string): path to raster containing the parameter
            pramn_1, the minimum aboveground ratio with zero biomass
        pramn_2_path (string): path to raster containing the parameter
            pramn_2, the minimum aboveground ratio with biomass greater
            than or equal to biomax
        aglivc_path (string): path to raster containing carbon in
            aboveground live biomass
        biomax_path (string): path to raster containing the parameter
            biomax, the biomass above which the ratio equals pramn_2
            or pramx_2
        pramx_1_path (string): path to raster containing the parameter
            pramx_1, the maximum aboveground ratio with zero biomass
        pramx_2_path (string): path to raster containing the parameter
            pramx_2, the maximum aboveground ratio with biomass greater
            than or equal to biomax
        prbmn_1_path (string): path to raster containing the parameter
            prbmn_1, intercept of regression to predict minimum
            belowground ratio from annual precipitation
        prbmn_2_path (string): path to raster containing the parameter
            prbmn_2, slope of regression to predict minimum belowground
            ratio from annual precipitation
        prbmx_1_path (string): path to raster containing the parameter
            prbmx_1, intercept of regression to predict maximum belowground
            ratio from annual precipitation
        prbmx_2_path (string): path to raster containing the parameter
            prbmx_2, slope of regression to predict maximum belowground
            ratio from annual precipitation
        annual_precip_path (string): path to annual precipitation raster
        month_reg (dict): map of key, path pairs giving paths to
            intermediate calculated values that are shared between
            submodels
        pft_i (int): plant functional type index
        iel (int): nutrient index (iel=1 indicates N, iel=2 indicates P)

    Modifies:
        rasters indicated by
            `month_reg['cercrp_min_above_<iel>_<pft_i>']`,
            `month_reg['cercrp_max_above_<iel>_<pft_i>']`,
            `month_reg['cercrp_min_below_<iel>_<pft_i>']`,
            `month_reg['cercrp_max_below_<iel>_<pft_i>']`,

    Returns:
        None
    """
    def calc_above_ratio(pra_1, pra_2, aglivc, biomax):
        """Calculate carbon to nutrient ratio for aboveground material.

        Parameters:
            pra_1 (numpy.ndarray): parameter, minimum or maximum ratio
                with zero biomass
            pra_2 (numpy.ndarray): parameter, minimum or maximum ratio
                with biomass greater than or equal to biomax
            aglivc (numpy.ndarray): state variable, carbon in aboveground
                live material
            biomax (numpy:ndarray): parameter, biomass above which the
                ratio equals pra_2

        Returns:
            cercrp_above, carbon to nutrient ratio for aboveground
                material
        """
        valid_mask = (
            (pra_1 != _IC_NODATA) &
            (pra_2 != _IC_NODATA) &
            (aglivc != aglivc_nodata) &
            (biomax != _IC_NODATA))

        cercrp_above = numpy.empty(pra_1.shape, dtype=numpy.float32)
        cercrp_above[:] = _TARGET_NODATA

        cercrp_above[valid_mask] = numpy.minimum(
            (pra_1[valid_mask] + (pra_2[valid_mask] - pra_1[valid_mask]) *
                2.5 * aglivc[valid_mask] / biomax[valid_mask]),
            pra_2[valid_mask])
        return cercrp_above

    def calc_below_ratio(prb_1, prb_2, annual_precip):
        """Calculate carbon to nutrient ratio for belowground material.

        Parameters:
            prb_1 (numpy.ndarray): parameter, intercept of regression
                to predict ratio from annual precipitation
            prb_2 (numpy.ndarray): parameter, slope of regression to
                predict ratio from annual precipitation
            annual_precip (numpy.ndarray): derived, precipitation in twelve
                months including the current month

        Returns:
            cercrp_below, carbon to nutrient ratio for belowground
                material
        """
        valid_mask = (
            (prb_1 != _IC_NODATA) &
            (prb_2 != _IC_NODATA) &
            (annual_precip != annual_precip_nodata))

        cercrp_below = numpy.empty(prb_1.shape, dtype=numpy.float32)
        cercrp_below[:] = _TARGET_NODATA

        cercrp_below[valid_mask] = (
            prb_1[valid_mask] +
            (prb_2[valid_mask] * annual_precip[valid_mask]))
        return cercrp_below

    aglivc_nodata = pygeoprocessing.get_raster_info(
        aglivc_path)['nodata'][0]
    annual_precip_nodata = pygeoprocessing.get_raster_info(
        annual_precip_path)['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            pramn_1_path, pramn_2_path, aglivc_path, biomax_path]],
        calc_above_ratio,
        month_reg['cercrp_min_above_{}_{}'.format(iel, pft_i)],
        gdal.GDT_Float32, _TARGET_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            pramx_1_path, pramx_2_path, aglivc_path, biomax_path]],
        calc_above_ratio,
        month_reg['cercrp_max_above_{}_{}'.format(iel, pft_i)],
        gdal.GDT_Float32, _TARGET_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            prbmn_1_path, prbmn_2_path, annual_precip_path]],
        calc_below_ratio,
        month_reg['cercrp_min_below_{}_{}'.format(iel, pft_i)],
        gdal.GDT_Float32, _TARGET_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            prbmx_1_path, prbmx_2_path, annual_precip_path]],
        calc_below_ratio,
        month_reg['cercrp_max_below_{}_{}'.format(iel, pft_i)],
        gdal.GDT_Float32, _TARGET_NODATA)


def calc_revised_fracrc(
        frtcindx_path, fracrc_p_path, totale_1_path, totale_2_path,
        demand_1_path, demand_2_path, h2ogef_1_path, cfrtcw_1_path,
        cfrtcw_2_path, cfrtcn_1_path, cfrtcn_2_path, fracrc_r_path):
    """
    Calculate revised fraction of carbon allocated to roots.

    The revised fraction of carbon allocated to roots includes the
    impacts of water and nutrient limitation. The method of the
    revised calculation depends on whether the plant functional
    type is modeled as a perennial plant or with the "Great Plains"
    equation of Parton et al. 1987, "Analysis of factors controlling soil
    organic matter levels in Great Plains grasslands", Soil Science
    Society of America Journal. Lines 96-104, cropDynC.f, froota.f

    Parameters:
        frtcindx_path (string): path to raster containing the parameter
            frtcindx
        fracrc_p_path (string): path to raster containing provisional
            fraction of carbon allocated to roots
        totale_1_path (string): path to raster containing total available
            nitrogen
        totale_2_path (string): path to raster containing total available
            phosphorus
        demand_1_path (string): path to raster containing nitrogen demand
        demand_2_path (string): path to raster containing phosphorus demand
        h2ogef_1_path (string): path to raster containing the limiting
            effect of water availability on growth
        cfrtcw_1_path (string): path to raster containing the parameter
            cfrtcw_1
        cfrtcw_2_path (string): path to raster containing the parameter
            cfrtcw_2
        cfrtcn_1_path (string): path to raster containing the parameter
            cfrtcn_1
        cfrtcn_2_path (string): path to raster containing the parameter
            cfrtcn_2
        fracrc_r_path (string): path to raster that should contain the
            result, revised fraction of carbon allocated to roots

    Modifies:
        the raster indicated by `fracrc_r_path`

    Returns:
        None
    """
    def calc_a2drat(totale, demand):
        """Calculate the ratio of available nutrient to nutrient demand.

        The ratio of nutrient available to demand for the nutrient is
        restricted to be between 0 and 1.

        Parameters:
            totale (numpy.ndarray): derived, nutrient available
            demand (numpy.ndarray): derived, demand for the nutrient

        Returns:
            a2drat, the ratio of available nutrient to demand, restricted
                to be between 0 and 1
        """
        valid_mask = (
            (totale != _TARGET_NODATA) &
            (demand != _TARGET_NODATA))

        a2drat = numpy.empty(totale.shape, dtype=numpy.float32)
        a2drat[:] = _TARGET_NODATA
        a2drat[valid_mask] = numpy.clip(
            totale[valid_mask] / demand[valid_mask], 0., 1.)
        return a2drat

    def calc_perennial_fracrc(
            h2ogef, cfrtcw_1, cfrtcw_2, a2drat_1, a2drat_2, cfrtcn_1,
            cfrtcn_2):
        """Calculate fraction C allocated to roots for a perennial plant.

        The fraction of carbon allocated to roots is determined by
        water availability, described by h2ogef, and nutrient availability,
        described by a2drat_1 for nitrogen and a2drat_2 for phosphorus.
        Lines 114-125 froota.f

        Parameters:
            h2ogef (numpy.ndarray): derived, the limiting factor of water
                availability on growth
            cfrtcw_1 (numpy.ndarray): parameter, the maximum fraction of
                carbon allocated to roots with maximum water stress
            cfrtcw_2 (numpy.ndarray): parameter, the minimum fraction of
                carbon allocated to roots with no water stress
            a2drat_1 (numpy.ndarray): derived, the ratio of available
                nitrogen to nitrogen demand, restricted to be between 0
                and 1
            a2drat_2 (numpy.ndarray): derived, the ratio of available
                phosphorus to phosphorus demand, restricted to be between
                0 and 1
            cfrtcn_1 (numpy.ndarray): parameter, maximum fraction of
                carbon allocated to roots with maximum nutrient stress
            cfrtcn_2 (numpy.ndarray): parameter, minimum fraction of
                carbon allocated to roots with no nutrient stress

        Returns:
            fracrc_perennial, revised fraction of C allocated to roots for
                a perennial plant
        """
        valid_mask = (
            (h2ogef != _TARGET_NODATA) &
            (cfrtcw_1 != _IC_NODATA) &
            (cfrtcw_2 != _IC_NODATA) &
            (a2drat_1 != _TARGET_NODATA) &
            (a2drat_2 != _TARGET_NODATA) &
            (cfrtcn_1 != _IC_NODATA) &
            (cfrtcn_2 != _IC_NODATA))

        h2oeff = numpy.empty(h2ogef.shape, dtype=numpy.float32)
        h2oeff[:] = _TARGET_NODATA
        h2oeff[valid_mask] = (
            (cfrtcw_2[valid_mask] - cfrtcw_1[valid_mask]) *
            (h2ogef[valid_mask] - 1.) + cfrtcw_2[valid_mask])

        ntreff_1 = numpy.empty(h2ogef.shape, dtype=numpy.float32)
        ntreff_1[:] = _TARGET_NODATA
        ntreff_1[valid_mask] = (
            (cfrtcn_2[valid_mask] - cfrtcn_1[valid_mask]) *
            (a2drat_1[valid_mask] - 1.0) + cfrtcn_2[valid_mask])

        ntreff_2 = numpy.empty(h2ogef.shape, dtype=numpy.float32)
        ntreff_2[:] = _TARGET_NODATA
        ntreff_1[valid_mask] = (
            (cfrtcn_2[valid_mask] - cfrtcn_1[valid_mask]) *
            (a2drat_2[valid_mask] - 1.0) + cfrtcn_2[valid_mask])

        ntreff = numpy.empty(h2ogef.shape, dtype=numpy.float32)
        ntreff[:] = _TARGET_NODATA
        ntreff[valid_mask] = numpy.maximum(
            ntreff_1[valid_mask], ntreff_2[valid_mask])

        fracrc_perennial = numpy.empty(
            h2ogef.shape, dtype=numpy.float32)
        fracrc_perennial[:] = _TARGET_NODATA
        fracrc_perennial[valid_mask] = numpy.minimum(
            numpy.maximum(h2oeff[valid_mask], ntreff[valid_mask]), 0.99)
        return fracrc_perennial

    def revised_fracrc_op(frtcindx, fracrc_p, fracrc_perennial):
        """Calculate revised fraction of carbon allocated to roots.

        The revised fraction of carbon allocated to roots is calculated
        according to the parameter frtcindx.  If frtcindx=0 (use the "Great
        Plains equation"), the revised fraction is equal to the provisional
        fraction.  If frtcindx=1 (a perennial plant), the revised fraction
        is calculated from water and nutrient stress.

        Parameters:
            frtcindx (numpy.ndarray): parameter, indicates whether revised
                fraction of carbon allocated to roots should follow the
                "Great Plains equation" or the algorithm for a perennial
                plant
            fracrc_p (numpy.ndarray): derived, provisional fraction of
                carbon allocated to roots
            fracrc_perennial (numpy.ndarray): derived, fraction of
                carbon allocated to roots for a perennial plant

        Returns:
            fracrc_r, revised fraction of carbon allocated to roots
        """
        valid_mask = (
            (frtcindx != _IC_NODATA) &
            (fracrc_p != _TARGET_NODATA) &
            (fracrc_perennial != _TARGET_NODATA))

        fracrc_r = numpy.empty(frtcindx.shape, dtype=numpy.float32)
        fracrc_r[:] = _TARGET_NODATA
        fracrc_r[valid_mask] = numpy.where(
            frtcindx[valid_mask] == 0, fracrc_p[valid_mask],
            fracrc_perennial[valid_mask])
        return fracrc_r

    # temporary intermediate rasters for calculating revised fracrc
    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {}
    for val in ['a2drat_1', 'a2drat_2', 'fracrc_perennial']:
        temp_val_dict[val] = os.path.join(
            temp_dir, '{}.tif'.format(val))

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [totale_1_path, demand_1_path]],
        calc_a2drat, temp_val_dict['a2drat_1'], gdal.GDT_Float32,
        _TARGET_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [totale_2_path, demand_2_path]],
        calc_a2drat, temp_val_dict['a2drat_2'], gdal.GDT_Float32,
        _TARGET_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            h2ogef_1_path, cfrtcw_1_path, cfrtcw_2_path,
            temp_val_dict['a2drat_1'], temp_val_dict['a2drat_2'],
            cfrtcn_1_path, cfrtcn_2_path]],
        calc_perennial_fracrc, temp_val_dict['fracrc_perennial'],
        gdal.GDT_Float32, _TARGET_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            frtcindx_path, fracrc_p_path,
            temp_val_dict['fracrc_perennial']]],
        revised_fracrc_op, fracrc_r_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def grazing_effect_on_aboveground_production(
        tgprod, fracrc, flgrem, grzeff):
    """Adjust aboveground production with the impact of grazing.

    Removal of biomass by herbivores directly impacts potential
    aboveground production according to the amount of biomass removed
    and the parameter grzeff, which acts as a switch to determine the
    effect. If grzeff=0, 3, or 4, aboveground production is not
    changed. If grzeff=1 or 6, production decreases linearly with
    biomass removed; if grzeff=2 or 5, biomass removed has a quadratic
    impact on production. Grazrst.f

    Parameters:
        tgprod (numpy.ndarray): derived, total potential biomass
            production restricted by water and nutrient availability
        fracrc (numpy.ndarray): derived, fraction of carbon allocated
            to roots according to water and nutrient availability
        flgrem (numpy.ndarray): derived, fraction of live biomass
            removed by grazing in previous monthly step
        grzeff (numpy.ndarray): parameter, the effect of defoliation on
            production and root:shoot ratio

    Returns:
        agprod, aboveground production impacted by grazing
    """
    valid_mask = (
        (tgprod != _TARGET_NODATA) &
        (fracrc != _TARGET_NODATA) &
        (flgrem != _IC_NODATA) &
        (grzeff != _IC_NODATA))

    agprod_prior = numpy.empty(tgprod.shape, dtype=numpy.float32)
    agprod_prior[:] = _TARGET_NODATA
    agprod_prior[valid_mask] = (
        tgprod[valid_mask] * (1. - fracrc[valid_mask]))

    linear_effect = numpy.empty(tgprod.shape, dtype=numpy.float32)
    linear_effect[:] = _TARGET_NODATA
    linear_effect[valid_mask] = numpy.maximum(
        (1. - (2.21*flgrem[valid_mask])) * agprod_prior[valid_mask],
        0.02)

    quadratic_effect = numpy.empty(tgprod.shape, dtype=numpy.float32)
    quadratic_effect[:] = _TARGET_NODATA
    quadratic_effect[valid_mask] = numpy.maximum(
        (1. + 2.6*flgrem[valid_mask] -
            (5.83*(numpy.power(flgrem[valid_mask], 2)))) *
        agprod_prior[valid_mask],
        0.02)

    no_effect_mask = (valid_mask & numpy.isin(grzeff, [0, 3, 4]))
    linear_mask = (valid_mask & numpy.isin(grzeff, [1, 6]))
    quadratic_mask = (valid_mask & numpy.isin(grzeff, [2, 5]))

    agprod = numpy.empty(tgprod.shape, dtype=numpy.float32)
    agprod[:] = _TARGET_NODATA
    agprod[no_effect_mask] = agprod_prior[no_effect_mask]
    agprod[linear_mask] = linear_effect[linear_mask]
    agprod[quadratic_mask] = quadratic_effect[quadratic_mask]
    return agprod


def grazing_effect_on_root_shoot(fracrc, flgrem, grzeff, gremb):
    """Adjust root:shoot ratio according to the impact of grazing.

    Removal of biomass by herbivores directly impacts the root:shoot
    ratio of production according to the amount of biomass removed and
    the parameter grzeff, which acts as a switch to determine the
    effect. If grzeff=0 or 1, the root:shoot ratio is not changed.
    If grzeff=2 or 3, biomass removed has a quadratic impact on the
    root:shoot ratio. If grzeff=4, 5, or 6, biomass removed has a
    linear effect on the root:shoot ratio. The parameter gremb
    multiplies the linear impact of grazing when grzeff=4, 5 or 6.
    Grzrst.f

    Parameters:
        fracrc (numpy.ndarray): derived, fraction of carbon allocated
            to roots according to water and nutrient availability
        flgrem (numpy.ndarray): derived, fraction of live biomass
            removed by grazing in previous monthly step
        grzeff (numpy.ndarray): parameter, the effect of defoliation on
            production and root:shoot ratio
        grzemb (numpy.ndarray): parameter, grazing effect multiplier

    Returns:
        rtsh, root:shoot ratio impacted by grazing
    """
    valid_mask = (
        (fracrc != _TARGET_NODATA) &
        (flgrem != _TARGET_NODATA) &
        (grzeff != _IC_NODATA) &
        (gremb != _IC_NODATA))

    rtsh_prior = numpy.empty(fracrc.shape, dtype=numpy.float32)
    rtsh_prior[:] = _TARGET_NODATA
    rtsh_prior[valid_mask] = (
        fracrc[valid_mask] / (1. - fracrc[valid_mask]))

    quadratic_effect = numpy.empty(fracrc.shape, dtype=numpy.float32)
    quadratic_effect[:] = _TARGET_NODATA
    quadratic_effect[valid_mask] = (
        rtsh_prior[valid_mask] + 3.05 * flgrem[valid_mask] -
        11.78 * numpy.power(flgrem[valid_mask], 2))

    linear_effect = numpy.empty(fracrc.shape, dtype=numpy.float32)
    linear_effect[:] = _TARGET_NODATA
    linear_effect[valid_mask] = (
        1. - (flgrem[valid_mask] * gremb[valid_mask]))

    no_effect_mask = (valid_mask & numpy.isin(grzeff, [0, 1]))
    quadratic_mask = (valid_mask & numpy.isin(grzeff, [2, 3]))
    linear_mask = (valid_mask & numpy.isin(grzeff, [4, 5, 6]))

    rtsh = numpy.empty(fracrc.shape, dtype=numpy.float32)
    rtsh[:] = _TARGET_NODATA
    rtsh[no_effect_mask] = rtsh_prior[no_effect_mask]
    rtsh[quadratic_mask] = quadratic_effect[quadratic_mask]
    rtsh[linear_mask] = linear_effect[linear_mask]
    return rtsh


def calc_tgprod_final(rtsh, agprod):
    """Calculate final total potential production.

    Final total potential production is calculated from aboveground
    production impacted by grazing and the final root:shoot ratio
    impacted by grazing.

    Parameters:
        rtsh (numpy.ndarray): derived, final root:shoot ratio impacted
            by grazing
        agprod (numpy.ndarray): derived, final aboveground potential
            production impacted by grazing

    Returns:
        tgprod, final total potential production
    """
    valid_mask = (
        (rtsh != _TARGET_NODATA) &
        (agprod != _TARGET_NODATA))
    tgprod = numpy.empty(rtsh.shape, dtype=numpy.float32)
    tgprod[:] = _TARGET_NODATA
    tgprod[valid_mask] = (
        agprod[valid_mask] + (rtsh[valid_mask] * agprod[valid_mask]))
    return tgprod


def calc_final_tgprod_rtsh(
        tgprod_pot_prod_path, fracrc_path, flgrem_path, grzeff_path,
        gremb_path, tgprod_path, rtsh_path):
    """Calculate final potential production and root:shoot ratio.

    Final potential production and root:shoot ratio include the impact of
    grazing. First calculate final aboveground production including the
    impact of grazing; then calculate rtsh, the final root:shoot ratio
    including the impact of grazing; then calculate tgprod, final total
    potential production, from final aboveground production and final
    root:shoot ratio. Grazrst.f

    Parameters:
        tgprod_pot_prod_path (string): path to raster containing total
            potential biomass production restricted by water and nutrient
            availability, prior to effects of grazing
        fracrc_path (string): path to raster containing the fraction of
            carbon production allocated to roots according to restriction
            by water and nutrient availability, prior to effects of
            grazing
        flgrem_path (string): path to raster containing the fraction of
            live aboveground biomass removed by herbivores according to
            diet selection in the previous step
        grzeff_path (string): path to raster containing the parameter
            grzeff, the effect of defolation on production and root:shoot
            ratio
        gremb_path (string): path to raster containing the parameter
            gremb, the grazing effect multiplier
        tgprod_path (string): path to raster containing final total
            potential production
        rtsh_path (string): path to raster containing final root:shoot
            ratio

    Modifies:
        The raster indicated by tgprod_path
        The raster indicated by rtsh_path

    Returns:
        None
    """
    # temporary intermediate rasters for grazing effect
    temp_dir = tempfile.mkdtemp()
    agprod_path = os.path.join(temp_dir, 'agprod.tif')

    # grazing effect on aboveground production
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in
            tgprod_pot_prod_path, fracrc_path, flgrem_path,
            grzeff_path],
        grazing_effect_on_aboveground_production,
        agprod_path, gdal.GDT_Float32, _TARGET_NODATA)
    # grazing effect on final root:shoot ratio
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in
            fracrc_path, flgrem_path, grzeff_path, gremb_path],
        grazing_effect_on_root_shoot, rtsh_path,
        gdal.GDT_Float32, _TARGET_NODATA)
    # final total potential production
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in rtsh_path, agprod_path],
        calc_tgprod_final, tgprod_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _root_shoot_ratio(
        aligned_inputs, site_param_table, current_month, pft_id_set,
        veg_trait_table, sv_reg, year_reg, month_reg):
    """Calculate final potential production and root:shoot ratio.

    Final potential production and root:shoot ratio is calculated according
    to nutrient availability and demand for the nutrient, and the impact of
    defoliation by herbivores. CropDynC.f

    Parameters:
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including the site spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        current_month (int): month of the year, such that current_month=1
            indicates January
        pft_id_set (set): set of integers identifying plant functional types
        veg_trait_table (dict): map of pft id to dictionaries containing
            plant functional type parameters
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month
        year_reg (dict): map of key, path pairs giving paths to rasters that
            are modified once per year, including annual precipitation
        month_reg (dict): map of key, path pairs giving paths to intermediate
            calculated values that are shared between submodels

    Modifies:
        The raster indicated by `month_reg['tgprod_<PFT>']` for each
            plant functional type (PFT)
        The raster indicated by `month_reg['rtsh_<PFT>']` for each
            plant functional type (PFT)

    Returns:
        None
    """
    # if growth does not occur this month for all PFTs,
    # skip the rest of the function
    do_PFT = []
    for pft_index in pft_id_set:
        if str(current_month) in veg_trait_table[pft_index]['growth_months']:
            do_PFT.append(pft_index)
        else:
            month_reg['tgprod_{}_path'.format(pft_index)] = None
            month_reg['rtsh_{}_path'.format(pft_index)] = None
    if not do_PFT:
        return

    # temporary intermediate rasters for root:shoot submodel
    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {}
    for pft_i in do_PFT:
        for val in ['fracrc_p', 'fracrc']:
            temp_val_dict['{}_{}'.format(val, pft_i)] = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
        for iel in [1, 2]:
            for val in ['eavail', 'demand']:
                temp_val_dict[
                    '{}_{}_{}'.format(val, pft_i, iel)] = os.path.join(
                        temp_dir, '{}_{}_{}.tif'.format(val, pft_i, iel))

    # temporary parameter rasters for root:shoot submodel
    param_val_dict = {}
    # site-level parameters
    for val in [
            'bgppa', 'bgppb', 'agppa', 'agppb', 'favail_1', 'favail_4',
            'favail_5', 'favail_6']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (aligned_inputs['site_index'], 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)
    # PFT-level parameters
    for pft_i in do_PFT:
        for val in [
                'frtcindx', 'cfrtcw_1', 'cfrtcw_2', 'cfrtcn_1', 'cfrtcn_2',
                'biomax', 'cfrtcw_1', 'cfrtcw_2', 'cfrtcn_1', 'cfrtcn_2',
                'grzeff', 'gremb']:
            target_path = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
            param_val_dict['{}_{}'.format(val, pft_i)] = target_path
            fill_val = veg_trait_table[pft_i][val]
            pygeoprocessing.new_raster_from_base(
                aligned_inputs['site_index'], target_path, gdal.GDT_Float32,
                [_IC_NODATA], fill_value_list=[fill_val])
        for iel in [1, 2]:
            for val in [
                    'pramn_1', 'pramn_2', 'pramx_1', 'pramx_2', 'prbmn_1',
                    'prbmn_2', 'prbmx_1', 'prbmx_2']:
                target_path = os.path.join(
                    temp_dir, '{}_{}_{}.tif'.format(val, pft_i, iel))
                param_val_dict[
                    '{}_{}_{}'.format(val, pft_i, iel)] = target_path
                fill_val = veg_trait_table[pft_i]['{}_{}'.format(val, iel)]
                pygeoprocessing.new_raster_from_base(
                    aligned_inputs['site_index'], target_path,
                    gdal.GDT_Float32, [_IC_NODATA], fill_value_list=[fill_val])

    # the parameter favail_2 must be calculated from current mineral N in
    # surface layer
    param_val_dict['favail_2'] = os.path.join(temp_dir, 'favail_2.tif')
    _calc_favail_P(sv_reg, param_val_dict)
    for pft_i in pft_id_set:
        # fracrc_p, provisional fraction of C allocated to roots
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in
                year_reg['annual_precip_path'],
                param_val_dict['frtcindx_{}'.format(pft_i)],
                param_val_dict['bgppa'],
                param_val_dict['bgppb'],
                param_val_dict['agppa'],
                param_val_dict['agppb'],
                param_val_dict['cfrtcw_1_{}'.format(pft_i)],
                param_val_dict['cfrtcw_2_{}'.format(pft_i)],
                param_val_dict['cfrtcn_1_{}'.format(pft_i)],
                param_val_dict['cfrtcn_2_{}'.format(pft_i)]],
            calc_provisional_fracrc,
            temp_val_dict['fracrc_p_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)
        for iel in [1, 2]:
            # persistent ratios used here and in plant growth submodel
            calc_ce_ratios(
                param_val_dict['pramn_1_{}_{}'.format(pft_i, iel)],
                param_val_dict['pramn_2_{}_{}'.format(pft_i, iel)],
                sv_reg['aglivc_{}_path'.format(pft_i)],
                param_val_dict['biomax_{}'.format(pft_i)],
                param_val_dict['pramx_1_{}_{}'.format(pft_i, iel)],
                param_val_dict['pramx_2_{}_{}'.format(pft_i, iel)],
                param_val_dict['prbmn_1_{}_{}'.format(pft_i, iel)],
                param_val_dict['prbmn_2_{}_{}'.format(pft_i, iel)],
                param_val_dict['prbmx_1_{}_{}'.format(pft_i, iel)],
                param_val_dict['prbmx_2_{}_{}'.format(pft_i, iel)],
                year_reg['annual_precip_path'], month_reg,
                pft_i, iel)
            # eavail_iel, available nutrient
            _calc_available_nutrient(
                pft_i, iel, veg_trait_table[pft_i], sv_reg, site_param_table,
                aligned_inputs['site_index'],
                param_val_dict['favail_{}'.format(iel)],
                month_reg['tgprod_pot_prod_{}'.format(pft_i)],
                temp_val_dict['eavail_{}_{}'.format(pft_i, iel)])
            # demand_iel, demand for the nutrient
            _calc_nutrient_demand(
                month_reg['tgprod_pot_prod_{}'.format(pft_i)],
                temp_val_dict['fracrc_p_{}'.format(pft_i)],
                month_reg['cercrp_min_above_{}_{}'.format(iel, pft_i)],
                month_reg['cercrp_min_below_{}_{}'.format(iel, pft_i)],
                temp_val_dict['demand_{}_{}'.format(pft_i, iel)])
        # revised fraction of carbon allocated to roots
        calc_revised_fracrc(
            param_val_dict['frtcindx_{}'.format(pft_i)],
            temp_val_dict['fracrc_p_{}'.format(pft_i)],
            temp_val_dict['eavail_{}_1'.format(pft_i)],
            temp_val_dict['eavail_{}_2'.format(pft_i)],
            temp_val_dict['demand_{}_1'.format(pft_i)],
            temp_val_dict['demand_{}_2'.format(pft_i)],
            month_reg['h2ogef_1_{}'.format(pft_i)],
            param_val_dict['cfrtcw_1_{}'.format(pft_i)],
            param_val_dict['cfrtcw_2_{}'.format(pft_i)],
            param_val_dict['cfrtcn_1_{}'.format(pft_i)],
            param_val_dict['cfrtcn_2_{}'.format(pft_i)],
            temp_val_dict['fracrc_{}'.format(pft_i)])
        # final potential production and root:shoot ratio accounting for
        # impacts of grazing
        # for now: TODO how to store flgrem and fdgrem?
        flgrem_path = os.path.join(temp_dir, 'flgrem.tif')
        pygeoprocessing.new_raster_from_base(
            param_val_dict['grzeff_{}'.format(pft_i)], flgrem_path,
            gdal.GDT_Float32, [_TARGET_NODATA], fill_value_list=[0.05])
        calc_final_tgprod_rtsh(
            month_reg['tgprod_pot_prod_{}'.format(pft_i)],
            temp_val_dict['fracrc_{}'.format(pft_i)],
            flgrem_path,
            param_val_dict['grzeff_{}'.format(pft_i)],
            param_val_dict['gremb_{}'.format(pft_i)],
            month_reg['tgprod_{}'.format(pft_i)],
            month_reg['rtsh_{}'.format(pft_i)])

    # clean up temporary files
    shutil.rmtree(temp_dir)
