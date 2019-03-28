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

# temporary directory to store intermediate files
PROCESSING_DIR = None

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
    'strlig_1_path': 'strlig_1.tif',
    'strlig_2_path': 'strlig_2.tif',
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
    'parent_2_path': 'parent_2.tif',
    'occlud_path': 'occlud.tif',
    'som1e_1_2_path': 'som1e_1_2.tif',
    'som1e_2_2_path': 'som1e_2_2.tif',
    'som2e_1_2_path': 'som2e_1_2.tif',
    'som2e_2_2_path': 'som2e_2_2.tif',
    'som3e_2_path': 'som3e_2.tif',
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
    'avh2o_3_path': 'avh2o_3.tif',
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
    'snow_path': 'snow.tif',
    'snlq_path': 'snlq.tif',
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

# site-level values that are updated once per year
_YEARLY_FILES = {
    'annual_precip_path': 'annual_precip.tif',
    'baseNdep_path': 'baseNdep.tif',
    }

# pft-level values that are updated once per year
_YEARLY_PFT_FILES = ['pltlig_above', 'pltlig_below']

# intermediate values for each plant functional type that are shared
# between submodels, but do not need to be saved as output
_PFT_INTERMEDIATE_VALUES = [
    'h2ogef_1', 'tgprod_pot_prod',
    'cercrp_min_above_1', 'cercrp_min_above_2',
    'cercrp_max_above_1', 'cercrp_max_above_2',
    'cercrp_min_below_1', 'cercrp_min_below_2',
    'cercrp_max_below_1', 'cercrp_max_below_2',
    'tgprod', 'rtsh']

# intermediate site-level values that are shared between submodels,
# but do not need to be saved as output
_SITE_INTERMEDIATE_VALUES = ['amov_2', 'snowmelt', 'bgwfunc']

# Target nodata is for general rasters that are positive, and _IC_NODATA are
# for rasters that are any range
_TARGET_NODATA = -1.0
_IC_NODATA = float(numpy.finfo('float32').min)
# SV_NODATA is for state variables
_SV_NODATA = -1.0


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
            (args['site_param_spatial_index_path'], 1)):
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
    for pft_i in pft_id_set:
        for sv in _PFT_STATE_VARIABLES:
            pft_sv_dict['{}_{}_path'.format(
                sv, pft_i)] = '{}_{}.tif'.format(sv, pft_i)

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

    # temporary directory for intermediate files
    PROCESSING_DIR = os.path.join(args['workspace_dir'], "temporary_files")
    if not os.path.exists(PROCESSING_DIR):
        os.makedirs(PROCESSING_DIR)

    # set up a dictionary that uses the same keys as
    # 'base_align_raster_path_id_map' to point to the clipped/resampled
    # rasters to be used in raster calculations for the model.
    aligned_raster_dir = os.path.join(
        args['workspace_dir'], 'aligned_inputs')
    os.makedirs(aligned_raster_dir)
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
        # set _SV_NODATA from initial rasters
        state_var_nodata = set([])
        # align initial state variables to resampled inputs
        resample_initial_path_map = {}
        for sv in _SITE_STATE_VARIABLE_FILES.iterkeys():
            sv_path = os.path.join(
                args['initial_conditions_dir'],
                _SITE_STATE_VARIABLE_FILES[sv])
            state_var_nodata.update(
                set([pygeoprocessing.get_raster_info(sv_path)['nodata'][0]]))
            resample_initial_path_map[sv] = sv_path
            if not os.path.exists(sv_path):
                missing_initial_values.append(sv_path)
        for pft_i in pft_id_set:
            for sv in _PFT_STATE_VARIABLES:
                sv_key = '{}_{}_path'.format(sv, pft_i)
                sv_path = os.path.join(
                    args['initial_conditions_dir'],
                    '{}_{}.tif'.format(sv, pft_i))
                state_var_nodata.update(
                    set([pygeoprocessing.get_raster_info(sv_path)['nodata']
                        [0]]))
                resample_initial_path_map[sv_key] = sv_path
                if not os.path.exists(sv_path):
                    missing_initial_values.append(sv_path)
        if missing_initial_values:
            raise ValueError(
                "Couldn't find the following required initial values: " +
                "\n\t".join(missing_initial_values))
        if len(state_var_nodata) > 1:
            raise ValueError(
                "Initial state variable rasters contain >1 nodata value")
        _SV_NODATA = list(state_var_nodata)[0]

        # align and resample initialization rasters with inputs
        sv_dir = os.path.join(args['workspace_dir'], 'state_variables_m-1')
        os.makedirs(sv_dir)
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
            ['near'] * len(initial_path_list),
            target_pixel_size, 'intersection',
            base_vector_path_list=[args['aoi_path']])
        sv_reg = aligned_initial_path_map
    else:
        LOGGER.info("initial conditions not supplied")
        LOGGER.info("aligning base raster inputs")
        pygeoprocessing.align_and_resize_raster_stack(
            source_input_path_list, aligned_input_path_list,
            ['near'] * len(aligned_inputs),
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
    year_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    year_reg = dict(
        [(key, os.path.join(year_dir, path)) for key, path in
            _YEARLY_FILES.iteritems()])
    for pft_i in pft_id_set:
        for file in _YEARLY_PFT_FILES:
            year_reg['{}_{}'.format(file, pft_i)] = os.path.join(
                year_dir, '{}_{}.tif'.format(file, pft_i))

    # make monthly directory for monthly intermediate parameters that are
    # shared between submodels, but do not need to be saved as output
    month_temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    month_reg = {}
    for pft_i in pft_id_set:
        for val in _PFT_INTERMEDIATE_VALUES:
            month_reg['{}_{}'.format(
                val, pft_i)] = os.path.join(
                month_temp_dir, '{}_{}.tif'.format(val, pft_i))
    for val in _SITE_INTERMEDIATE_VALUES:
        month_reg[val] = os.path.join(month_temp_dir, '{}.tif'.format(val))

    # Main simulation loop
    # for each step in the simulation
    for month_index in xrange(n_months):
        if (month_index % 12) == 0:
            # Update yearly quantities
            _yearly_tasks(
                aligned_inputs, site_param_table, veg_trait_table, month_index,
                year_reg, pft_id_set)

        current_month = (starting_month + month_index - 1) % 12 + 1
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

        _potential_production(
            aligned_inputs, site_param_table, current_month, month_index,
            pft_id_set, veg_trait_table, prev_sv_reg, pp_reg, month_reg)

        _root_shoot_ratio(
            aligned_inputs, site_param_table, current_month, pft_id_set,
            veg_trait_table, prev_sv_reg, year_reg, month_reg)

        _soil_water(
            aligned_inputs, site_param_table, veg_trait_table, current_month,
            month_index, prev_sv_reg, sv_reg, pp_reg, month_reg, pft_id_set)

        _decomposition(
            aligned_inputs, current_month, month_index, site_param_table,
            year_reg, month_reg, prev_sv_reg, sv_reg, pp_reg)

        _death_and_partition(
            'stded', aligned_inputs, site_param_table, current_month,
            prev_sv_reg, sv_reg, year_reg, pft_id_set, veg_trait_table)

        _death_and_partition(
            'bgliv', aligned_inputs, site_param_table, current_month,
            prev_sv_reg, sv_reg, year_reg, pft_id_set, veg_trait_table)

        _shoot_senescence(
            pft_id_set, veg_trait_table, prev_sv_reg, sv_reg, month_reg,
            current_month)

        _new_growth(
            pft_id_set, aligned_inputs, site_param_table, veg_trait_table,
            sv_reg, month_reg, current_month)


def raster_multiplication(
        raster1, raster1_nodata, raster2, raster2_nodata, target_path,
        target_path_nodata):
    """Multiply raster1 by raster2.

    Multiply raster1 by raster2 element-wise. In any pixel where raster1 or
    raster2 is nodata, the result is nodata. The result is always of float
    datatype.

    Modifies:
        the raster indicated by `target_path`

    Returns:
        None
    """
    def raster_multiply_op(raster1, raster2):
        """Multiply two rasters."""
        valid_mask = (
            (~numpy.isclose(raster1, raster1_nodata)) &
            (~numpy.isclose(raster2, raster2_nodata)))
        result = numpy.empty(raster1.shape, dtype=numpy.float32)
        result[:] = target_path_nodata
        result[valid_mask] = raster1[valid_mask] * raster2[valid_mask]
        return result
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [raster1, raster2]],
        raster_multiply_op, target_path, gdal.GDT_Float32,
        target_path_nodata)


def raster_division(
        raster1, raster1_nodata, raster2, raster2_nodata, target_path,
        target_path_nodata):
    """Divide raster1 by raster2.

    Divide raster1 by raster2 element-wise. In any pixel where raster1 or
    raster2 is nodata, the result is nodata. The result is always of float
    datatype.

    Modifies:
        the raster indicated by `target_path`

    Returns:
        None
    """
    def raster_divide_op(raster1, raster2):
        """Divide raster1 by raster2."""
        valid_mask = (
            (~numpy.isclose(raster1, raster1_nodata)) &
            (~numpy.isclose(raster2, raster2_nodata)))
        raster1 = raster1.astype(numpy.float32)
        raster2 = raster2.astype(numpy.float32)

        result = numpy.empty(raster1.shape, dtype=numpy.float32)
        result[:] = target_path_nodata
        error_mask = ((raster1 != 0) & (raster2 == 0.) & valid_mask)
        zero_mask = ((raster1 == 0.) & (raster2 == 0.) & valid_mask)
        nonzero_mask = ((raster2 != 0.) & valid_mask)
        result[error_mask] = target_path_nodata
        result[zero_mask] = 0.
        result[nonzero_mask] = raster1[nonzero_mask] / raster2[nonzero_mask]
        return result
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [raster1, raster2]],
        raster_divide_op, target_path, gdal.GDT_Float32,
        target_path_nodata)


def raster_list_sum(
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
            numpy.isclose(numpy.array(raster_list), input_nodata), axis=0)
        sum_of_rasters = numpy.sum(raster_list, axis=0)
        sum_of_rasters[invalid_mask] = target_nodata
        return sum_of_rasters

    def raster_sum_op_nodata_remove(*raster_list):
        """Add the rasters in raster_list, treating nodata as zero."""
        for r in raster_list:
            numpy.place(r, numpy.isclose(r, input_nodata), [0])
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


def raster_sum(
        raster1, raster1_nodata, raster2, raster2_nodata, target_path,
        target_nodata, nodata_remove=False):
    """Add raster 1 and raster2.

    Add raster1 and raster2, allowing nodata values in the rasters to
    propagate to the result or treating nodata as zero.

    Parameters:
        raster1 (string): path to one raster operand
        raster1_nodata (float or int): nodata value in raster1
        raster2 (string): path to second raster operand
        raster2_nodata (float or int): nodata value in raster2
        target_path (string): path to location to store the sum
        target_nodata (float or int): nodata value for the result raster
        nodata_remove (bool): if true, treat nodata values in input
            rasters as zero. If false, the sum in a pixel where any
            input raster is nodata is nodata.

    Modifies:
        the raster indicated by `target_path`

    Returns:
        None
    """
    def raster_sum_op(raster1, raster2):
        """Add raster1 and raster2 without removing nodata values."""
        valid_mask = (
            (~numpy.isclose(raster1, raster1_nodata)) &
            (~numpy.isclose(raster2, raster2_nodata)))
        result = numpy.empty(raster1.shape, dtype=numpy.float32)
        result[:] = target_nodata
        result[valid_mask] = raster1[valid_mask] + raster2[valid_mask]
        return result

    def raster_sum_op_nodata_remove(raster1, raster2):
        """Subtract raster2 from raster1, treating nodata as zero."""
        numpy.place(raster1, numpy.isclose(raster1, raster1_nodata), [0])
        numpy.place(raster2, numpy.isclose(raster2, raster2_nodata), [0])
        result = raster1 + raster2
        return result

    if nodata_remove:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [raster1, raster2]],
            raster_sum_op_nodata_remove, target_path, gdal.GDT_Float32,
            target_nodata)
    else:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [raster1, raster2]],
            raster_sum_op, target_path, gdal.GDT_Float32,
            target_nodata)


def raster_difference(
        raster1, raster1_nodata, raster2, raster2_nodata, target_path,
        target_nodata, nodata_remove=False):
    """Subtract raster2 from raster1.

    Subtract raster2 from raster1 element-wise, allowing nodata values in the
    rasters to propagate to the result or treating nodata as zero.

    Parameters:
        raster1 (string): path to raster from which to subtract raster2
        raster1_nodata (float or int): nodata value in raster1
        raster2 (string): path to raster which should be subtracted from
            raster1
        raster2_nodata (float or int): nodata value in raster2
        target_path (string): path to location to store the difference
        target_nodata (float or int): nodata value for the result raster
        nodata_remove (bool): if true, treat nodata values in input
            rasters as zero. If false, the difference in a pixel where any
            input raster is nodata is nodata.

    Modifies:
        the raster indicated by `target_path`

    Returns:
        None
    """
    def raster_difference_op(raster1, raster2):
        """Subtract raster2 from raster1 without removing nodata values."""
        valid_mask = (
            (~numpy.isclose(raster1, raster1_nodata)) &
            (~numpy.isclose(raster2, raster2_nodata)))
        result = numpy.empty(raster1.shape, dtype=numpy.float32)
        result[:] = target_nodata
        result[valid_mask] = raster1[valid_mask] - raster2[valid_mask]
        return result

    def raster_difference_op_nodata_remove(raster1, raster2):
        """Subtract raster2 from raster1, treating nodata as zero."""
        numpy.place(raster1, numpy.isclose(raster1, raster1_nodata), [0])
        numpy.place(raster2, numpy.isclose(raster2, raster2_nodata), [0])
        result = raster1 - raster2
        return result

    if nodata_remove:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [raster1, raster2]],
            raster_difference_op_nodata_remove, target_path, gdal.GDT_Float32,
            target_nodata)
    else:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [raster1, raster2]],
            raster_difference_op, target_path, gdal.GDT_Float32,
            target_nodata)


def reclassify_nodata(target_path, new_nodata_value):
    """Reclassify the nodata value of a raster to a new value.

    Convert all areas of nodata in the target raster to the new nodata
    value, which must be an integer.

    Parameters:
        target_path (string): path to target raster
        new_nodata_value (integer): new value to set as nodata

    Modifies:
        the raster indicated by `target_path`

    Returns:
        None
    """
    def reclassify_op(target_raster):
        reclassified_raster = numpy.copy(target_raster)
        reclassify_mask = (target_raster == previous_nodata_value)
        reclassified_raster[reclassify_mask] = new_nodata_value
        return reclassified_raster

    fd, temp_path = tempfile.mkstemp(dir=PROCESSING_DIR)
    shutil.copyfile(target_path, temp_path)
    previous_nodata_value = pygeoprocessing.get_raster_info(
        target_path)['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(temp_path, 1)], reclassify_op, target_path, gdal.GDT_Float32,
        new_nodata_value)

    # clean up
    os.close(fd)
    os.remove(temp_path)


def weighted_state_variable_sum(
        sv, sv_reg, aligned_inputs, pft_id_set, weighted_sum_path):
    """Calculate weighted sum of state variable across plant functional types.

    To sum a state variable across PFTs within a grid cell, the state variable
    must be weighted by the percent cover of each PFT inside the grid cell.
    First multiply the state variable by its percent cover, and then add up the
    weighted products.

    Parameters:
        sv (string): state variable to be summed across plant functional types
        sv_reg (dict): map of key, path pairs giving paths to state variables,
            including sv, the state variable to be summed
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including percent cover of each plant
            functional type
        pft_id_set (set): set of integers identifying plant functional types
        weighted_sum_path (string): path to raster that should contain the
            weighted sum across PFTs

    Modifies:
        the raster indicated by `weighted_sum_path`

    Returns:
        None
    """
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for pft_i in pft_id_set:
        val = '{}_weighted'.format(sv)
        temp_val_dict['{}_{}'.format(val, pft_i)] = os.path.join(
            temp_dir, '{}_{}.tif'.format(val, pft_i))

    weighted_path_list = []
    for pft_i in pft_id_set:
        target_path = temp_val_dict['{}_weighted_{}'.format(sv, pft_i)]
        pft_nodata = pygeoprocessing.get_raster_info(
            aligned_inputs['pft_{}'.format(pft_i)])['nodata'][0]
        raster_multiplication(
            sv_reg['{}_{}_path'.format(sv, pft_i)], _SV_NODATA,
            aligned_inputs['pft_{}'.format(pft_i)], pft_nodata,
            target_path, _TARGET_NODATA)
        weighted_path_list.append(target_path)
    raster_list_sum(
        weighted_path_list, _TARGET_NODATA, weighted_sum_path, _TARGET_NODATA,
        nodata_remove=True)

    # clean up temporary files
    shutil.rmtree(temp_dir)


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
            (~numpy.isclose(som1c_2, _SV_NODATA)) &
            (~numpy.isclose(som2c_2, _SV_NODATA)) &
            (~numpy.isclose(som3c, _SV_NODATA)) &
            (~numpy.isclose(bulkd, bulkd_nodata)) &
            (edepth != _IC_NODATA))
        ompc[valid_mask] = (
            (som1c_2[valid_mask] + som2c_2[valid_mask] +
                som3c[valid_mask]) * 1.724 /
            (10000. * bulkd[valid_mask] * edepth[valid_mask]))
        return ompc

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
            (~numpy.isclose(sand, sand_nodata)) &
            (~numpy.isclose(silt, silt_nodata)) &
            (~numpy.isclose(clay, clay_nodata)) &
            (ompc != _TARGET_NODATA) &
            (~numpy.isclose(bulkd, bulkd_nodata)))
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
            (~numpy.isclose(sand, sand_nodata)) &
            (~numpy.isclose(silt, silt_nodata)) &
            (~numpy.isclose(clay, clay_nodata)) &
            (ompc != _TARGET_NODATA) &
            (~numpy.isclose(bulkd, bulkd_nodata)))
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
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
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
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
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
            (~numpy.isclose(sand, sand_nodata)))
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
            (~numpy.isclose(sand, sand_nodata)))
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
            (~numpy.isclose(clay, clay_nodata)))
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
            (~numpy.isclose(clay, clay_nodata)))
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
            (~numpy.isclose(sand, sand_nodata)))
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


def _aboveground_ratio(anps, tca, pcemic_1, pcemic_2, pcemic_3):
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

    Returns:
        agdrat, the C/<iel> ratio of new material
    """
    valid_mask = (
        (~numpy.isclose(anps, _SV_NODATA)) &
        (~numpy.isclose(tca, _SV_NODATA)) &
        (pcemic_1 != _IC_NODATA) &
        (pcemic_2 != _IC_NODATA) &
        (pcemic_3 != _IC_NODATA))

    cemicb = numpy.empty(anps.shape, dtype=numpy.float32)
    cemicb[:] = _IC_NODATA
    cemicb[valid_mask] = (
        (pcemic_2[valid_mask] - pcemic_1[valid_mask]) /
        pcemic_3[valid_mask])

    econt = numpy.empty(anps.shape, dtype=numpy.float32)
    econt[:] = _TARGET_NODATA
    econt[valid_mask] = 0

    decompose_mask = ((tca > 0.) & valid_mask)
    econt[decompose_mask] = anps[decompose_mask] / (tca[decompose_mask] * 2.5)

    agdrat = numpy.empty(anps.shape, dtype=numpy.float32)
    agdrat[:] = _TARGET_NODATA
    agdrat[valid_mask] = pcemic_2[valid_mask]

    compute_mask = ((econt <= pcemic_3) & valid_mask)
    agdrat[compute_mask] = (
        pcemic_1[compute_mask] + econt[compute_mask] * cemicb[compute_mask])
    return agdrat


def _belowground_ratio(aminrl, varat_1_iel, varat_2_iel, varat_3_iel):
    """Calculate C/<iel> ratios of decomposing belowground material.

    This ratio is used to test whether there is sufficient <iel> (N or P)
    in soil metabolic material to decompose.  Bgdrat.f

    Parameters:
        aminrl (numpy.ndarray): derived, average surface mineral <iel>
        varat_1_iel (numpy.ndarray): parameter, maximum C/<iel> ratio for
            newly decomposed material
        varat_2_iel (numpy.ndarray): parameter, minimum C/<iel> ratio
        varat_3_iel (numpy.ndarray): parameter, amount of <iel> present
            when minimum ratio applies

    Returns:
        bgdrat, the C/<iel> ratio of new material
    """
    valid_mask = (
        (~numpy.isclose(aminrl, _SV_NODATA)) &
        (varat_1_iel != _IC_NODATA) &
        (varat_2_iel != _IC_NODATA) &
        (varat_3_iel != _IC_NODATA))

    bgdrat = numpy.empty(aminrl.shape, dtype=numpy.float32)
    bgdrat[:] = _TARGET_NODATA
    bgdrat[valid_mask] = (
        (1. - aminrl[valid_mask] / varat_3_iel[valid_mask]) *
        (varat_1_iel[valid_mask] - varat_2_iel[valid_mask]) +
        varat_2_iel[valid_mask])

    max_mask = ((aminrl <= 0) & valid_mask)
    bgdrat[max_mask] = varat_1_iel[max_mask]

    min_mask = ((aminrl > varat_3_iel) & valid_mask)
    bgdrat[min_mask] = varat_2_iel[min_mask]
    return bgdrat


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
    # temporary parameter rasters for structural ratios calculations
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
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
            (~numpy.isclose(struce_1, _SV_NODATA)) &
            (~numpy.isclose(strucc_1, _SV_NODATA)) &
            (rad1p_1 != _IC_NODATA) &
            (rad1p_2 != _IC_NODATA) &
            (rad1p_3 != _IC_NODATA) &
            (pcemic1_2 != _IC_NODATA) &
            (rnewas1 != _TARGET_NODATA))

        rnewas2 = _aboveground_ratio(
            struce_1, strucc_1, pcemic2_1, pcemic2_2, pcemic2_3)

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
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sv_reg['struce_1_{}_path'.format(iel)],
                sv_reg['strucc_1_path'],
                param_val_dict['pcemic1_1_{}'.format(iel)],
                param_val_dict['pcemic1_2_{}'.format(iel)],
                param_val_dict['pcemic1_3_{}'.format(iel)]]],
            _aboveground_ratio, pp_reg['rnewas_{}_1_path'.format(iel)],
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
                pp_reg['rnewas_{}_1_path'.format(iel)]]],
            calc_rnewas_som2, pp_reg['rnewas_{}_2_path'.format(iel)],
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
        aligned_inputs, site_param_table, veg_trait_table, month_index,
        year_reg, pft_id_set):
    """Calculate quantities that remain static for 12 months.

    These quantities are annual precipitation, annual atmospheric N
    deposition, and the fraction of plant residue which is lignin for each pft.
    Century also calculates non-symbiotic soil N fixation once yearly, but here
    those were moved to monthly tasks. Century uses precipitation in the future
    12 months (prcgrw) to predict root:shoot ratios, but here we instead use
    the sum of monthly precipitation in 12 months including the current one, if
    data for 12 future months are not available.
    Lines 79-82, 164 Eachyr.f

    Parameters:
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including monthly precipitation and site
            spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        veg_trait_table (dict): map of pft id to dictionaries containing
            plant functional type parameters
        month_index (int): current monthly step, relative to 0 so that
            month_index=0 at first monthly time step
        year_reg (dict): map of key, path pairs giving paths to the annual
            precipitation and N deposition rasters
        pft_id_set (set): set of integers identifying plant functional types

    Modifies:
        the rasters indicated by:
            year_reg['annual_precip_path']
            year_reg['baseNdep_path']
            year_reg['pltlig_above_<pft>'] for each pft
            year_reg['pltlig_below_<pft>'] for each pft

    Returns:
        None

    Raises:
        ValueError if fewer than 12 monthly precipitation rasters can be found
    """
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

    def calc_pltlig(fligni_1_lyr, fligni_2_lyr, prcann):
        """Calculate the fraction of residue that is lignin. Cmplig.f

        This fraction is used to calculate the fraction of residue (i.e.,
        incoming litter from fall of standing dead or incoming soil from death
        of roots) that is partitioned to metabolic vs structural pools. It is
        calculated once per year from annual precipitation and fixed
        parameters.

        Parameters:
            fligni_1_lyr (numpy.ndarray): parameter, intercept for regression
                predicting lignin content fraction from rainfall
            fligni_2_lyr (numpy.ndarray): parameter, slope for regression
                predicting lignin content fraction from rainfall
            prcann (numpy.ndarray): derived, annual precipitation

        Returns:
            pltlig_lyr, fraction of residue that is lignin
        """
        valid_mask = (
            (fligni_1_lyr != _IC_NODATA) &
            (fligni_2_lyr != _IC_NODATA) &
            (prcann != _TARGET_NODATA))

        pltlig = numpy.empty(fligni_1_lyr.shape, dtype=numpy.float32)
        pltlig[:] = _TARGET_NODATA
        pltlig[valid_mask] = (
            fligni_1_lyr[valid_mask] + fligni_2_lyr[valid_mask] *
            prcann[valid_mask])
        pltlig[valid_mask] = numpy.clip(pltlig[valid_mask], 0.02, 0.5)
        return pltlig

    offset = -12
    annual_precip_rasters = []
    while len(annual_precip_rasters) < 12:
        offset += 1
        if offset == 12:
            raise ValueError("Insufficient precipitation rasters were found")
        precip_month = month_index + offset
        try:
            annual_precip_rasters.append(
                    aligned_inputs['precip_%d' % precip_month])
        except KeyError:
            continue

    precip_nodata = set([])
    for precip_raster in annual_precip_rasters:
        precip_nodata.update(
            set([pygeoprocessing.get_raster_info(precip_raster)['nodata'][0]]))
    if len(precip_nodata) > 1:
        raise ValueError("Precipitation rasters include >1 nodata value")
    precip_nodata = list(precip_nodata)[0]

    raster_list_sum(
        annual_precip_rasters, precip_nodata, year_reg['annual_precip_path'],
        _TARGET_NODATA)

    # intermediate parameter rasters for this operation
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    param_val_dict = {}
    for val in['epnfa_1', 'epnfa_2']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (aligned_inputs['site_index'], 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)
    for val in ['fligni_1_1', 'fligni_2_1', 'fligni_1_2', 'fligni_2_2']:
        for pft_i in pft_id_set:
            target_path = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
            param_val_dict['{}_{}'.format(val, pft_i)] = target_path
            fill_val = veg_trait_table[pft_i][val]
            pygeoprocessing.new_raster_from_base(
                aligned_inputs['site_index'], target_path, gdal.GDT_Float32,
                [_IC_NODATA], fill_value_list=[fill_val])

    # calculate base N deposition
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['epnfa_1'], param_val_dict['epnfa_2'],
            year_reg['annual_precip_path']]],
        calc_base_N_dep, year_reg['baseNdep_path'], gdal.GDT_Float32,
        _TARGET_NODATA)

    for pft_i in pft_id_set:
        # fraction of surface residue that is lignin
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                param_val_dict['fligni_1_1_{}'.format(pft_i)],
                param_val_dict['fligni_2_1_{}'.format(pft_i)],
                year_reg['annual_precip_path']]],
            calc_pltlig, year_reg['pltlig_above_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

        # fraction of soil residue that is lignin
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                param_val_dict['fligni_1_2_{}'.format(pft_i)],
                param_val_dict['fligni_2_2_{}'.format(pft_i)],
                year_reg['annual_precip_path']]],
            calc_pltlig, year_reg['pltlig_below_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

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
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    latitude_raster_path = os.path.join(temp_dir, 'latitude.tif')
    pygeoprocessing.new_raster_from_base(
        template_raster, latitude_raster_path, gdal.GDT_Float32,
        [_IC_NODATA])
    latitude_raster = gdal.OpenEx(latitude_raster_path, gdal.GA_Update)
    target_band = latitude_raster.GetRasterBand(1)
    base_raster_info = pygeoprocessing.get_raster_info(template_raster)
    geotransform = base_raster_info['geotransform']
    for offset_map, raster_block in pygeoprocessing.iterblocks(
            (template_raster, 1)):
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
            (~numpy.isclose(max_temp, maxtmp_nodata)) &
            (~numpy.isclose(min_temp, mintmp_nodata)) &
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
        pft_id_set, veg_trait_table, prev_sv_reg, pp_reg, month_reg):
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
        prev_sv_reg (dict): map of key, path pairs giving paths to state
            variables for the previous month
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
    for pft_i in pft_id_set:
        if str(current_month) in veg_trait_table[pft_i]['growth_months']:
            do_PFT.append(pft_i)
    if not do_PFT:
        return

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
            (~numpy.isclose(maxtmp, maxtmp_nodata)) &
            (pmxtmp != _IC_NODATA) &
            (~numpy.isclose(mintmp, mintmp_nodata)) &
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
            avh2o_1 (numpy.ndarray): state variable, water available to this
                plant functional type for growth
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
            (~numpy.isclose(avh2o_1, _SV_NODATA)) &
            (~numpy.isclose(precip, precip_nodata)) &
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

        intcpt = (
            pprpts_1[valid_mask] + (pprpts_2[valid_mask] * wc[valid_mask]))
        slope = 1. / (pprpts_3[valid_mask] - intcpt)

        h2ogef_1 = numpy.empty(pevap.shape, dtype=numpy.float32)
        h2ogef_1[:] = _TARGET_NODATA
        h2ogef_1[valid_mask] = (
            1.0 + slope *
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
            (~numpy.isclose(strucc_1, _SV_NODATA)) &
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
                production (g biomass)
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
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
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

    maxtmp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['max_temp_{}'.format(current_month)])['nodata'][0]
    mintmp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['min_temp_{}'.format(current_month)])['nodata'][0]
    precip_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['precip_{}'.format(month_index)])['nodata'][0]

    # calculate intermediate quantities that do not differ between PFTs:
    # sum of aglivc (standing live biomass) and stdedc (standing dead biomass)
    # across PFTs, weighted by % cover of each PFT
    for sv in ['aglivc', 'stdedc']:
        weighted_sum_path = temp_val_dict['sum_{}'.format(sv)]
        weighted_state_variable_sum(
            sv, prev_sv_reg, aligned_inputs, pft_id_set, weighted_sum_path)

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
        # potprd, the limiting effect of temperature
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['ctemp'],
                param_val_dict['ppdf_1_{}'.format(pft_i)],
                param_val_dict['ppdf_2_{}'.format(pft_i)],
                param_val_dict['ppdf_3_{}'.format(pft_i)],
                param_val_dict['ppdf_4_{}'.format(pft_i)]]],
            calc_potprd, temp_val_dict['potprd_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

        # h2ogef_1, the limiting effect of soil water availability
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['pevap'],
                prev_sv_reg['avh2o_1_{}_path'.format(pft_i)],
                aligned_inputs['precip_{}'.format(month_index)],
                pp_reg['wc_path'],
                param_val_dict['pprpts_1'],
                param_val_dict['pprpts_2'],
                param_val_dict['pprpts_3']]],
            calc_h2ogef_1, month_reg['h2ogef_1_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

        # biof, the limiting effect of obstruction
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['sum_stdedc'],
                temp_val_dict['sum_aglivc'],
                prev_sv_reg['strucc_1_path'],
                param_val_dict['pmxbio'],
                param_val_dict['biok5_{}'.format(pft_i)]]],
            calc_biof, temp_val_dict['biof_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)

        # total potential production
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                param_val_dict['prdx_1_{}'.format(pft_i)],
                temp_val_dict['shwave'],
                temp_val_dict['potprd_{}'.format(pft_i)],
                month_reg['h2ogef_1_{}'.format(pft_i)],
                temp_val_dict['biof_{}'.format(pft_i)]]],
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
            (~numpy.isclose(minerl_1_1, _SV_NODATA)) &
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

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            sv_reg['minerl_1_1_path'],
            param_val_dict['favail_4'],
            param_val_dict['favail_5'],
            param_val_dict['favail_6']]],
        favail_P_op, param_val_dict['favail_2'],
        gdal.GDT_Float32, _IC_NODATA)


def _calc_avail_mineral_nutrient(pft_param_dict, sv_reg, iel, target_path):
    """Calculate one mineral nutrient available to one plant functional type.

    The mineral nutrient available to a plant functional type is calculated
    from the mineral nutrient content of soil layers accessible by that
    plant function type.

    Parameters:
        pft_param_dict (dict): map of key, value pairs giving the values of
            parameters for this plant functional type
            (i.e., veg_trait_table[pft_i] for this pft_i)
        sv_reg (dict): map of key, path pairs giving paths to state
            variables for the current month
        iel (int): integer index for current nutrient (1=N, 2=P)
        target_path (string): path to raster to contain available mineral
            nutrient for this plant functional type and nutrient

    Modifies:
        the raster indicated by `target_path`

    Returns:
        None
    """
    nlay = int(pft_param_dict['nlaypg'])
    mineral_raster_list = [
        sv_reg['minerl_{}_{}_path'.format(lyr, iel)] for lyr in xrange(
            1, nlay + 1)]
    raster_list_sum(
        mineral_raster_list, _SV_NODATA, target_path, _TARGET_NODATA,
        nodata_remove=True)


def _calc_available_nutrient(
        pft_i, iel, pft_param_dict, sv_reg, site_param_table, site_index_path,
        availm_path, favail_path, tgprod_path, eavail_path):
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
        availm_path (string): path to raster containing available mineral
            nutrient for the given plant functional type and nutrient
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        favail_path (string): path to raster containing the appropriate value
            of the parameter favail. For nitrogen, this parameter is supplied
            directly as user input, but for phosphorus, it must be calculated
            from other parameters.
        tgprod_path (string): path to raster containing total potential
            production (g biomass)
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
            (~numpy.isclose(bglivc, _SV_NODATA)) &
            (riint != _IC_NODATA) &
            (availm != _TARGET_NODATA) &
            (favail != _IC_NODATA) &
            (~numpy.isclose(crpstg, _SV_NODATA)))

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
                potential production (g biomass)

        Returns:
            eavail, total N available including N fixed by the plant

        """
        valid_mask = (
            (eavail_prior != _TARGET_NODATA) &
            (snfxmx != _IC_NODATA) &
            (tgprod != _TARGET_NODATA))

        maxNfix = numpy.empty(eavail_prior.shape, dtype=numpy.float32)
        maxNfix[:] = _TARGET_NODATA
        maxNfix[valid_mask] = snfxmx[valid_mask] * (tgprod[valid_mask] / 2.5)

        eavail = numpy.empty(eavail_prior.shape, dtype=numpy.float32)
        eavail[:] = _TARGET_NODATA
        eavail[valid_mask] = eavail_prior[valid_mask] + maxNfix[valid_mask]
        return eavail

    # temporary intermediate rasters for calculating available nutrient
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
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

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            param_val_dict['rictrl'],
            sv_reg['bglivc_{}_path'.format(pft_i)],
            param_val_dict['riint'],
            availm_path, favail_path,
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
        [(path, 1) for path in [
            biomass_production_path, fraction_allocated_to_roots_path,
            cercrp_min_above_path, cercrp_min_below_path]],
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
            pramn_<iel>_1, the minimum aboveground ratio with zero biomass
        pramn_2_path (string): path to raster containing the parameter
            pramn_<iel>_2, the minimum aboveground ratio with biomass greater
            than or equal to biomax
        aglivc_path (string): path to raster containing carbon in
            aboveground live biomass
        biomax_path (string): path to raster containing the parameter
            biomax, the biomass above which the ratio equals pramn_2
            or pramx_2
        pramx_1_path (string): path to raster containing the parameter
            pramx_<iel>_1, the maximum aboveground ratio with zero biomass
        pramx_2_path (string): path to raster containing the parameter
            pramx_<iel>_2, the maximum aboveground ratio with biomass greater
            than or equal to biomax
        prbmn_1_path (string): path to raster containing the parameter
            prbmn_<iel>_1, intercept of regression to predict minimum
            belowground ratio from annual precipitation
        prbmn_2_path (string): path to raster containing the parameter
            prbmn_<iel>_2, slope of regression to predict minimum belowground
            ratio from annual precipitation
        prbmx_1_path (string): path to raster containing the parameter
            prbmx_<iel>_1, intercept of regression to predict maximum
            belowground ratio from annual precipitation
        prbmx_2_path (string): path to raster containing the parameter
            prbmx_<iel>_2, slope of regression to predict maximum belowground
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
            (~numpy.isclose(aglivc, _SV_NODATA)) &
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
            (annual_precip != _TARGET_NODATA))

        cercrp_below = numpy.empty(prb_1.shape, dtype=numpy.float32)
        cercrp_below[:] = _TARGET_NODATA

        cercrp_below[valid_mask] = (
            prb_1[valid_mask] +
            (prb_2[valid_mask] * annual_precip[valid_mask]))
        return cercrp_below

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
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
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
    quadratic_effect[valid_mask] = (
        (1. + 2.6*flgrem[valid_mask] -
            (5.83*(numpy.power(flgrem[valid_mask], 2)))) *
        agprod_prior[valid_mask])
    quadratic_effect[valid_mask] = numpy.maximum(
        quadratic_effect[valid_mask], 0.02)

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
            potential production (g biomass)
        rtsh_path (string): path to raster containing final root:shoot
            ratio of potential production

    Modifies:
        The raster indicated by tgprod_path
        The raster indicated by rtsh_path

    Returns:
        None
    """
    # temporary intermediate rasters for grazing effect
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    agprod_path = os.path.join(temp_dir, 'agprod.tif')

    # grazing effect on aboveground production
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            tgprod_pot_prod_path, fracrc_path, flgrem_path,
            grzeff_path]],
        grazing_effect_on_aboveground_production,
        agprod_path, gdal.GDT_Float32, _TARGET_NODATA)
    # grazing effect on final root:shoot ratio
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            fracrc_path, flgrem_path, grzeff_path, gremb_path]],
        grazing_effect_on_root_shoot, rtsh_path,
        gdal.GDT_Float32, _TARGET_NODATA)
    # final total potential production
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [rtsh_path, agprod_path]],
        calc_tgprod_final, tgprod_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _root_shoot_ratio(
        aligned_inputs, site_param_table, current_month, pft_id_set,
        veg_trait_table, prev_sv_reg, year_reg, month_reg):
    """Calculate final potential production and root:shoot ratio.

    Final potential biomass production and root:shoot ratio is calculated
    according to nutrient availability and demand for the nutrient, and the
    impact of defoliation by herbivores. CropDynC.f

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
        prev_sv_reg (dict): map of key, path pairs giving paths to state
            variables for the previous month
        year_reg (dict): map of key, path pairs giving paths to rasters that
            are modified once per year, including annual precipitation
        month_reg (dict): map of key, path pairs giving paths to intermediate
            calculated values that are shared between submodels

    Modifies:
        The raster indicated by `month_reg['tgprod_<PFT>']`, total potential
            production (g biomass) for each plant functional type (PFT)
        The raster indicated by `month_reg['rtsh_<PFT>']` for each
            plant functional type (PFT)

    Returns:
        None

    """
    # if growth does not occur this month for all PFTs,
    # skip the rest of the function
    do_PFT = []
    for pft_i in pft_id_set:
        if str(current_month) in veg_trait_table[pft_i]['growth_months']:
            do_PFT.append(pft_i)
        else:
            month_reg['tgprod_{}'.format(pft_i)] = None
            month_reg['rtsh_{}'.format(pft_i)] = None
    if not do_PFT:
        return

    # temporary intermediate rasters for root:shoot submodel
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for pft_i in do_PFT:
        for val in ['fracrc_p', 'fracrc', 'availm']:
            temp_val_dict['{}_{}'.format(val, pft_i)] = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
        for iel in [1, 2]:
            for val in ['eavail', 'demand']:
                temp_val_dict[
                    '{}_{}_{}'.format(val, iel, pft_i)] = os.path.join(
                        temp_dir, '{}_{}_{}.tif'.format(val, iel, pft_i))

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
        for val in [
                'pramn_1_1', 'pramn_1_2', 'pramx_1_1', 'pramx_1_2',
                'prbmn_1_1', 'prbmn_1_2', 'prbmx_1_1', 'prbmx_1_2',
                'pramn_2_1', 'pramn_2_2', 'pramx_2_1', 'pramx_2_2',
                'prbmn_2_1', 'prbmn_2_2', 'prbmx_2_1', 'prbmx_2_2']:
            target_path = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
            param_val_dict[
                '{}_{}'.format(val, pft_i)] = target_path
            fill_val = veg_trait_table[pft_i][val]
            pygeoprocessing.new_raster_from_base(
                aligned_inputs['site_index'], target_path,
                gdal.GDT_Float32, [_IC_NODATA], fill_value_list=[fill_val])

    # the parameter favail_2 must be calculated from current mineral N in
    # surface layer
    param_val_dict['favail_2'] = os.path.join(temp_dir, 'favail_2.tif')
    _calc_favail_P(prev_sv_reg, param_val_dict)
    for pft_i in pft_id_set:
        # fracrc_p, provisional fraction of C allocated to roots
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                year_reg['annual_precip_path'],
                param_val_dict['frtcindx_{}'.format(pft_i)],
                param_val_dict['bgppa'],
                param_val_dict['bgppb'],
                param_val_dict['agppa'],
                param_val_dict['agppb'],
                param_val_dict['cfrtcw_1_{}'.format(pft_i)],
                param_val_dict['cfrtcw_2_{}'.format(pft_i)],
                param_val_dict['cfrtcn_1_{}'.format(pft_i)],
                param_val_dict['cfrtcn_2_{}'.format(pft_i)]]],
            calc_provisional_fracrc,
            temp_val_dict['fracrc_p_{}'.format(pft_i)],
            gdal.GDT_Float32, _TARGET_NODATA)
        for iel in [1, 2]:
            # persistent ratios used here and in plant growth submodel
            calc_ce_ratios(
                param_val_dict['pramn_{}_1_{}'.format(iel, pft_i)],
                param_val_dict['pramn_{}_2_{}'.format(iel, pft_i)],
                prev_sv_reg['aglivc_{}_path'.format(pft_i)],
                param_val_dict['biomax_{}'.format(pft_i)],
                param_val_dict['pramx_{}_1_{}'.format(iel, pft_i)],
                param_val_dict['pramx_{}_2_{}'.format(iel, pft_i)],
                param_val_dict['prbmn_{}_1_{}'.format(iel, pft_i)],
                param_val_dict['prbmn_{}_2_{}'.format(iel, pft_i)],
                param_val_dict['prbmx_{}_1_{}'.format(iel, pft_i)],
                param_val_dict['prbmx_{}_2_{}'.format(iel, pft_i)],
                year_reg['annual_precip_path'], month_reg,
                pft_i, iel)
            # sum of mineral nutrient in accessible soil layers
            _calc_avail_mineral_nutrient(
                veg_trait_table[pft_i], prev_sv_reg, iel,
                temp_val_dict['availm_{}'.format(pft_i)])
            # eavail_iel, available nutrient
            _calc_available_nutrient(
                pft_i, iel, veg_trait_table[pft_i], prev_sv_reg,
                site_param_table, aligned_inputs['site_index'],
                temp_val_dict['availm_{}'.format(pft_i)],
                param_val_dict['favail_{}'.format(iel)],
                month_reg['tgprod_pot_prod_{}'.format(pft_i)],
                temp_val_dict['eavail_{}_{}'.format(iel, pft_i)])
            # demand_iel, demand for the nutrient
            _calc_nutrient_demand(
                month_reg['tgprod_pot_prod_{}'.format(pft_i)],
                temp_val_dict['fracrc_p_{}'.format(pft_i)],
                month_reg['cercrp_min_above_{}_{}'.format(iel, pft_i)],
                month_reg['cercrp_min_below_{}_{}'.format(iel, pft_i)],
                temp_val_dict['demand_{}_{}'.format(iel, pft_i)])
        # revised fraction of carbon allocated to roots
        calc_revised_fracrc(
            param_val_dict['frtcindx_{}'.format(pft_i)],
            temp_val_dict['fracrc_p_{}'.format(pft_i)],
            temp_val_dict['eavail_1_{}'.format(pft_i)],
            temp_val_dict['eavail_2_{}'.format(pft_i)],
            temp_val_dict['demand_1_{}'.format(pft_i)],
            temp_val_dict['demand_2_{}'.format(pft_i)],
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


def _snow(
        site_index_path, site_param_table, precip_path, tave_path,
        max_temp_path, min_temp_path, prev_snow_path, prev_snlq_path,
        current_month, snowmelt_path, snow_path, snlq_path,
        inputs_after_snow_path, pet_rem_path):
    """Account for precipitation as snow and snowmelt from snowpack.

    Determine whether precipitation falls as snow. Track the fate of
    new and existing snowpack including evaporation and melting. Track the
    the remaining snowpack and liquid in snow and potential
    evapotranspiration remaining after evaporation of snow.  Snowcent.f

    Parameters:
        site_index_path (string): path to site spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        precip_path (string): path to raster containing precipitation for the
            current month
        tave_path (string): path to raster containing average temperature for
            the current month
        max_temp_path (string): path to raster containing maximum temperature
            for the current month
        min_temp_path (string): path to raster containing minimum temperature
            for the current month
        prev_snow_path (string): path to raster containing current snowpack
        prev_snlq_path (string): path to raster containing current liquid in
            snow
        current_month (int): current month of the year, such that month=0
            indicates January
        snow_path (string): path to raster to contain modified snowpack
        snlq_path (string): path to raster to contain modified liquid in snow
        inputs_after_snow_path (string): path to raster containing water inputs
            to the system after accounting for snow
        pet_rem_path (string): path to raster containing potential
            evapotranspiration remaining after any evaporation of snow

    Modifies:
        the raster indicated by `snowmelt_path`
        the raster indicated by `snow_path`
        the raster indicated by `snlq_path`
        the raster indicated by `inputs_after_snow_path`
        the raster indicated by `pet_rem_path`

    Returns:
        None
    """
    def calc_snow_moisture(return_type):
        """Calculate change in snow, pet, snow liquid, and moisture inputs.

        Record changes in snowpack, liquid in snow, potential
        evapotranspiration energy, and liquid draining into soil from snow.

        Parameters:
            return_type (string): flag indicating whether modified snowpack,
                modified liquid in snow, modified potential evapotranspiration,
                or soil moisture inputs after snow should be returned

        Returns:
            the function `_calc_snow_moisture`
        """
        def _calc_snow_moisture(
                tave, precip, snow, snlq, pet, tmelt_1, tmelt_2, shwave):
            """Calculate the fate of moisture from snow.

            Calculate new snowfall or rain on snow. Calculate direct
            evaporation of snow and consumption of potential
            evapotranspiration energy. Calculate snowmelt and liquid draining
            from snow into the soil.

            Parameters:
                tave (numpy.ndarray): derived, average temperature
                precip (numpy.ndarray): input, precipitation for this month
                snow (numpy.ndarray): derived, existing snowpack prior to new
                    snowfall
                snlq (numpy.ndarray): derived, existing liquid in snowpack
                pet (numpy.ndarray): derived, potential evapotranspiration
                tmelt_1 (numpy.ndarray): parameter, minimum temperature above
                    which snow will melt
                tmelt_2 (numpy.ndarray): parameter, ratio between degrees above
                    the minimum temperature and cm of snow that will melt
                shwave (numpy.ndarray): derived, shortwave radiation outside
                    the atmosphere

            Returns:
                snowmelt if return_type is 'snowmelt'
                snow_revised if return_type is 'snow'
                snlq_revised if return_type is 'snlq'
                pet_revised if return_type is 'pet'
                inputs_after_snow if return_type is 'inputs_after_snow'
            """
            valid_mask = (
                (tave != _IC_NODATA) &
                (~numpy.isclose(precip, precip_nodata)) &
                (~numpy.isclose(snow, _SV_NODATA)) &
                (~numpy.isclose(snlq, _SV_NODATA)) &
                (pet != _TARGET_NODATA) &
                (tmelt_1 != _IC_NODATA) &
                (tmelt_2 != _IC_NODATA) &
                (shwave != _TARGET_NODATA))
            inputs_after_snow = numpy.empty(precip.shape, dtype=numpy.float32)
            inputs_after_snow[:] = _TARGET_NODATA
            inputs_after_snow[valid_mask] = precip[valid_mask]

            snowfall_mask = (valid_mask & (tave <= 0))
            snow[snowfall_mask] = (snow[snowfall_mask] + precip[snowfall_mask])
            inputs_after_snow[snowfall_mask] = 0.

            rain_on_snow_mask = (
                (valid_mask) &
                (tave > 0) &
                (snow > 0))
            snlq[rain_on_snow_mask] = (
                snlq[rain_on_snow_mask] + precip[rain_on_snow_mask])
            inputs_after_snow[rain_on_snow_mask] = 0.

            snowtot = numpy.zeros(snow.shape, dtype=numpy.float32)
            snowtot[valid_mask] = numpy.maximum(
                snow[valid_mask] + snlq[valid_mask], 0)

            evap_mask = (valid_mask & (snowtot > 0.))
            evsnow = numpy.zeros(snow.shape, dtype=numpy.float32)
            evsnow[evap_mask] = numpy.minimum(
                snowtot[evap_mask], pet[evap_mask] * 0.87)
            snow_revised = numpy.empty(snow.shape, dtype=numpy.float32)
            snow_revised[:] = _TARGET_NODATA
            snow_revised[valid_mask] = snow[valid_mask]
            snow_revised[evap_mask] = numpy.maximum(
                snow[evap_mask] - evsnow[evap_mask] *
                (snow[evap_mask] / snowtot[evap_mask]), 0.)
            snlq_revised = numpy.zeros(snow.shape, dtype=numpy.float32)
            snlq_revised[valid_mask] = snlq[valid_mask]
            snlq_revised[evap_mask] = numpy.maximum(
                snlq[evap_mask] - evsnow[evap_mask] *
                (snlq[evap_mask] / snowtot[evap_mask]), 0.)
            pet_revised = numpy.empty(snow.shape, dtype=numpy.float32)
            pet_revised[:] = _TARGET_NODATA
            pet_revised[valid_mask] = pet[valid_mask]
            pet_revised[evap_mask] = numpy.maximum(
                (pet[evap_mask] - evsnow[evap_mask] / 0.87), 0.)

            melt_mask = (valid_mask & (tave >= tmelt_1))
            snowmelt = numpy.zeros(snow.shape, dtype=numpy.float32)
            snowmelt[melt_mask] = numpy.clip(
                tmelt_2[melt_mask] * (tave[melt_mask] - tmelt_1[melt_mask]) *
                shwave[melt_mask], 0., snow_revised[melt_mask])
            snow_revised[melt_mask] = (
                snow_revised[melt_mask] - snowmelt[melt_mask])
            snlq_revised[melt_mask] = (
                snlq_revised[melt_mask] + snowmelt[melt_mask])

            drain_mask = (melt_mask & (snlq_revised > 0.5 * snow_revised))
            inputs_after_snow[drain_mask] = (
                snlq_revised[drain_mask] - 0.5 * snow_revised[drain_mask])
            snlq_revised[drain_mask] = (
                snlq_revised[drain_mask] - inputs_after_snow[drain_mask])

            if return_type == 'snowmelt':
                return snowmelt
            elif return_type == 'snow':
                return snow_revised
            elif return_type == 'snlq':
                return snlq_revised
            elif return_type == 'pet':
                return pet_revised
            else:
                return inputs_after_snow
        return _calc_snow_moisture

    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for val in ['shwave', 'pet']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))
    param_val_dict = {}
    for val in ['tmelt_1', 'tmelt_2', 'fwloss_4']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for (
                site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_val, target_path, gdal.GDT_Float32,
            _IC_NODATA)

    max_temp_nodata = pygeoprocessing.get_raster_info(
        max_temp_path)['nodata'][0]
    min_temp_nodata = pygeoprocessing.get_raster_info(
        min_temp_path)['nodata'][0]
    precip_nodata = pygeoprocessing.get_raster_info(
        precip_path)['nodata'][0]

    # solar radiation outside the atmosphere
    _shortwave_radiation(precip_path, current_month, temp_val_dict['shwave'])

    # pet, reference evapotranspiration modified by fwloss parameter
    _reference_evapotranspiration(
        max_temp_path, min_temp_path, temp_val_dict['shwave'],
        param_val_dict['fwloss_4'], temp_val_dict['pet'])

    # calculate snowmelt
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            tave_path, precip_path, prev_snow_path,
            prev_snlq_path, temp_val_dict['pet'],
            param_val_dict['tmelt_1'], param_val_dict['tmelt_2'],
            temp_val_dict['shwave']]],
        calc_snow_moisture('snowmelt'), snowmelt_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate change in snow
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            tave_path, precip_path, prev_snow_path,
            prev_snlq_path, temp_val_dict['pet'],
            param_val_dict['tmelt_1'], param_val_dict['tmelt_2'],
            temp_val_dict['shwave']]],
        calc_snow_moisture("snow"), snow_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate change in liquid in snow
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            tave_path, precip_path, prev_snow_path,
            prev_snlq_path, temp_val_dict['pet'],
            param_val_dict['tmelt_1'], param_val_dict['tmelt_2'],
            temp_val_dict['shwave']]],
        calc_snow_moisture("snlq"), snlq_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate change in potential evapotranspiration energy
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            tave_path, precip_path, prev_snow_path,
            prev_snlq_path, temp_val_dict['pet'],
            param_val_dict['tmelt_1'], param_val_dict['tmelt_2'],
            temp_val_dict['shwave']]],
        calc_snow_moisture("pet"), pet_rem_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate soil moisture inputs draining from snow after snowmelt
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            tave_path, precip_path, prev_snow_path,
            prev_snlq_path, temp_val_dict['pet'],
            param_val_dict['tmelt_1'], param_val_dict['tmelt_2'],
            temp_val_dict['shwave']]],
        calc_snow_moisture("inputs_after_snow"), inputs_after_snow_path,
        gdal.GDT_Float32, _TARGET_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _calc_aboveground_live_biomass(sum_aglivc, sum_tgprod):
    """Calculate aboveground live biomass for purposes of soil water.

    Live biomass impacts loss of moisture inputs through canopy
    interception and evapotranspiration. Because soil moisture is computed
    after potential production, but before actual growth of plants, some of
    the predicted growth in biomass (i.e., tgprod) is added here to
    existing standing live biomass (i.e., aglivc * 2.5; line 80,
    potprod.f, in Century).

    Parameters:
        sum_aglivc (numpy.ndarray): the sum of aglivc across plant
            functional types (pft), weighted by % cover of the pft
        sum_tgprod (numpy.ndarray): sum of tgprod, potential production
            limited by soil water, nutrient availability, and grazing,
            across pfts weighted by % cover of the pft

    Returns:
        aliv, aboveground live biomass for soil water submodel
    """
    valid_mask = (
        (sum_aglivc != _TARGET_NODATA) &
        (sum_tgprod != _TARGET_NODATA))
    aliv = numpy.empty(sum_aglivc.shape, dtype=numpy.float32)
    aliv[:] = _TARGET_NODATA
    aliv[valid_mask] = (
        sum_aglivc[valid_mask] * 2.5 + (0.25 * sum_tgprod[valid_mask]))
    return aliv


def _calc_standing_biomass(aliv, sum_stdedc):
    """Calculate total aboveground standing biomass for soil water.

    Total standing biomass impacts loss of moisture inputs by increasing
    total canopy interception and decreasing bare soil evaporation. It is
    the sum of live and dead standing biomass across plant functional
    types, bounded to be <= 800 g/m2.

    Parameters:
        aliv (numpy.ndarray): aboveground live biomass, calculated from
            aglivc and tgprod across plant functional types
        sum_stdedc (numpy.ndarray): aboveground standing dead C summed
            across plant functional types

    Returns:
        sd, total aboveground standing biomass for soil water.
    """
    valid_mask = (
        (aliv != _TARGET_NODATA) &
        (sum_stdedc != _TARGET_NODATA))
    sd = numpy.empty(aliv.shape, dtype=numpy.float32)
    sd[:] = _TARGET_NODATA
    sd[valid_mask] = numpy.minimum(
        aliv[valid_mask] + (sum_stdedc[valid_mask] * 2.5), 800.)
    return sd


def subtract_surface_losses(return_type):
    """Calculate surface losses to runoff and surface evaporation.

    Calculate the loss of surface moisture to runoff, canopy interception,
    and bare soil evaporation.

    Parameters:
        return_type (string): flag indicating whether soil moisture inputs
            after surface losses or total surface evaporation should be
            returned

    Returns:
        the function `_subtract_surface_losses`
    """
    def _subtract_surface_losses(
            inputs_after_snow, fracro, precro, snow, alit, sd, fwloss_1,
            fwloss_2, pet_rem):
        """Subtract moisture losses to runoff, interception, and evaporation.

        Of the surface water inputs from precipitation and snowmelt, some water
        is lost to runoff (line 113, H2olos.f). After runoff, some water is
        lost to canopy interception and bare soil evaporation, if there is no
        snow cover. Loss to canopy interception and bare soil evaporation is
        a function of live, standing dead, and surface litter biomass.  The
        total loss of moisture to interception and bare soil evaporation is
        bounded to be less than or equal to 40% of reference
        evapotranspiration.

        Parameters:
            inputs_after_snow (numpy.ndarray): derived, surface water inputs
                from precipitation and snowmelt, prior to runoff
            fracro (numpy.ndarray): parameter, fraction of surface water
                above precro that is lost to runoff
            precro (numpy.ndarray): parameter, amount of surface water that
                must be available for runoff to occur
            snow (numpy.ndarray): derived, current snowpack
            alit (numpy.ndarray): derived, biomass in surface litter
            sd (numpy.ndarray): derived, total standing biomass
            fwloss_1 (numpy.ndarray): parameter, scaling factor for
                interception and evaporation of precip by vegetation
            fwloss_2 (numpy.ndarray): parameter, scaling factor for bare soil
                evaporation of precip
            pet_rem (numpy.ndarray): derived, potential evaporation remaining
                after evaporation of snow

        Returns:
            inputs_after_surface, surface water inputs to soil after runoff
                and surface evaporation are subtracted, if return_type is
                'inputs_after_surface'
            absevap, bare soil evaporation, if return_type is 'absevap'
            evap_losses, total surface evaporation, if return_type is
                'evap_losses'
        """
        valid_mask = (
            (inputs_after_snow != _TARGET_NODATA) &
            (fracro != _IC_NODATA) &
            (precro != _IC_NODATA) &
            (snow != _TARGET_NODATA) &
            (alit != _TARGET_NODATA) &
            (sd != _TARGET_NODATA) &
            (fwloss_1 != _IC_NODATA) &
            (fwloss_2 != _IC_NODATA) &
            (pet_rem != _TARGET_NODATA))

        runoff = numpy.empty(inputs_after_snow.shape, dtype=numpy.float32)
        runoff[:] = _TARGET_NODATA
        runoff[valid_mask] = numpy.maximum(
            fracro[valid_mask] *
            (inputs_after_snow[valid_mask] - precro[valid_mask]), 0.)
        inputs_after_runoff = numpy.empty(
            inputs_after_snow.shape, dtype=numpy.float32)
        inputs_after_runoff[:] = _TARGET_NODATA
        inputs_after_runoff[valid_mask] = (
            inputs_after_snow[valid_mask] - runoff[valid_mask])

        evap_mask = (valid_mask & (snow <= 0))
        # loss to interception
        aint = numpy.zeros(inputs_after_snow.shape, dtype=numpy.float32)
        aint[evap_mask] = (
            (0.0003 * alit[evap_mask] + 0.0006 * sd[evap_mask]) *
            fwloss_1[evap_mask])
        # loss to bare soil evaporation
        absevap = numpy.empty(inputs_after_snow.shape, dtype=numpy.float32)
        absevap[:] = _TARGET_NODATA
        absevap[evap_mask] = (
            0.5 *
            numpy.exp((-0.002 * alit[evap_mask]) - (0.004 * sd[evap_mask])) *
            fwloss_2[evap_mask])
        # total losses to interception and evaporation
        evap_losses = numpy.empty(inputs_after_snow.shape, dtype=numpy.float32)
        evap_losses[:] = _TARGET_NODATA
        evap_losses[evap_mask] = (
            numpy.minimum(((absevap[evap_mask] + aint[evap_mask]) *
            inputs_after_runoff[evap_mask]), (0.4 * pet_rem[evap_mask])))
        # remaining inputs after evaporation
        inputs_after_surface = numpy.empty(
            inputs_after_snow.shape, dtype=numpy.float32)
        inputs_after_surface[:] = _TARGET_NODATA
        inputs_after_surface[valid_mask] = inputs_after_runoff[valid_mask]
        inputs_after_surface[evap_mask] = (
            inputs_after_runoff[evap_mask] - evap_losses[evap_mask])

        if return_type == 'inputs_after_surface':
            return inputs_after_surface
        elif return_type == 'absevap':
            return absevap
        elif return_type == 'evap_losses':
            return evap_losses
    return _subtract_surface_losses


def calc_potential_transpiration(return_type):
    """Calculate potential transpiration and evaporation from soil layer 1.

    Calculate potential transpiration (trap), potential evaporation from
    soil layer 1 (pevp), and initial transpiration water loss (tran).
    Remove the initial transpiration water loss from soil moisture inputs
    at this step.

    Parameters:
        return_type (string): flag indicating whether potential transpiration,
            potential evaporation from soil layer 1, or modified moisture
            inputs should be returned

    Returns:
        the function `_calc_potential_transpiration`
    """
    def _calc_potential_transpiration(
            pet_rem, evap_losses, tave, aliv, current_moisture_inputs):
        """Calculate potential water losses to transpiration.

        Calculate potential transpiration (trap), the total potential
        transpiration from all soil layers by plants. Calculate potential
        evaporation from soil layer 1 (pevp); this amount is calculated prior
        to transpiration but actually removed after water loss to transpiration
        from all soil layers has been accounted. Calculate actual transpiration
        (tran). Remove actual transpiration water losses from moisture inputs
        before distributing water to soil layers. This is necessary for a
        monthly time step to give plants in wet climates adequate access to
        water for transpiration.

        Parameters:
            pet_rem (numpy.ndarray): derived, potential evapotranspiration
                remaining after evaporation of snow
            evap_losses (numpy.ndarray): derived, total surface evaporation
            tave (numpy.ndarray): derived, average temperature
            aliv (numpy.ndarray): aboveground live biomass, calculated from
                aglivc and tgprod across plant functional types
            current_moisture_inputs (numpy.ndarray): derived, moisture inputs
                after surface losses

        Returns:
            trap if return_type is 'trap'
            pevp if return_type is 'pevp'
            modified_moisture_inputs if return_type is
                'modified_moisture_inputs'
        """
        valid_mask = (
            (pet_rem != _TARGET_NODATA) &
            (evap_losses != _TARGET_NODATA) &
            (tave != _TARGET_NODATA) &
            (aliv != _TARGET_NODATA) &
            (current_moisture_inputs != _TARGET_NODATA))
        trap = numpy.empty(pet_rem.shape, dtype=numpy.float32)
        trap[:] = _TARGET_NODATA
        trap[valid_mask] = pet_rem[valid_mask] - evap_losses[valid_mask]

        no_transpiration_mask = (valid_mask & (tave < 2))
        trap[no_transpiration_mask] = 0.

        transpiration_mask = (valid_mask & (tave >= 2))
        trap[transpiration_mask] = numpy.maximum(
            numpy.minimum(
                trap[transpiration_mask], pet_rem[transpiration_mask] *
                0.65 * (1 - numpy.exp(-0.02 * aliv[transpiration_mask]))), 0.)

        trap[valid_mask] = numpy.maximum(trap[valid_mask], 0.01)
        pevp = numpy.empty(pet_rem.shape, dtype=numpy.float32)
        pevp[:] = _TARGET_NODATA
        pevp[valid_mask] = numpy.maximum(
            pet_rem[valid_mask] - trap[valid_mask] - evap_losses[valid_mask],
            0.)

        tran = numpy.empty(pet_rem.shape, dtype=numpy.float32)
        tran[:] = _TARGET_NODATA
        tran[valid_mask] = numpy.minimum(
            trap[valid_mask] - 0.01, current_moisture_inputs[valid_mask])
        trap[valid_mask] = trap[valid_mask] - tran[valid_mask]

        modified_moisture_inputs = numpy.empty(
            pet_rem.shape, dtype=numpy.float32)
        modified_moisture_inputs[:] = _TARGET_NODATA
        modified_moisture_inputs[valid_mask] = (
            current_moisture_inputs[valid_mask] - tran[valid_mask])

        if return_type == 'trap':
            return trap
        elif return_type == 'pevp':
            return pevp
        elif return_type == 'modified_moisture_inputs':
            return modified_moisture_inputs
    return _calc_potential_transpiration


def distribute_water_to_soil_layer(return_type):
    """Distribute moisture inputs to one soil layer prior to transpiration.

    Soil moisture inputs after runoff, evaporation, and initial
    transpiration are distributed to soil layers sequentially according to the
    field capacity of the layer.  If moisture inputs exceed the field capacity
    of the layer, the remainder of moisture inputs move down to the next
    adjacent soil layer.

    Returns:
        the function `_distribute_water`
    """
    def _distribute_water(adep, afiel, asmos, current_moisture_inputs):
        """Revise soil moisture in this soil layer prior to transpiration.

        Moisture inputs coming into this soil layer are compared to the field
        capacity of the layer. If the field capacity is exceeded, the excess
        moisture moves from this layer to the next adjacent layer.

        Parameters:
            adep (numpy.ndarray): parameter, depth of this soil layer in cm
            afiel (numpy.ndarray): derived, field capacity of this layer
            asmos (numpy.ndarray): state variable, current soil moisture
                content of this soil layer
            current_moisture_inputs (numpy.ndarray): derived, moisture inputs
                added to this soil layer

        Returns:
            asmos_revised, revised soil moisture in this layer, if return_type
                is 'asmos_revised'
            amov, moisture flowing from this layer into the next, if
                return_type is 'amov'
        """
        valid_mask = (
            (adep != _IC_NODATA) &
            (afiel != _TARGET_NODATA) &
            (~numpy.isclose(asmos, _SV_NODATA)) &
            (current_moisture_inputs != _TARGET_NODATA))

        afl = numpy.empty(adep.shape, dtype=numpy.float32)
        afl[:] = _TARGET_NODATA
        afl[valid_mask] = adep[valid_mask] * afiel[valid_mask]

        asmos_interm = numpy.empty(adep.shape, dtype=numpy.float32)
        asmos_interm[:] = _TARGET_NODATA
        asmos_interm[valid_mask] = (
            asmos[valid_mask] + current_moisture_inputs[valid_mask])

        amov = numpy.empty(adep.shape, dtype=numpy.float32)
        amov[:] = _TARGET_NODATA

        exceeded_mask = (valid_mask & (asmos_interm > afl))
        amov[exceeded_mask] = asmos_interm[exceeded_mask]

        asmos_revised = numpy.empty(adep.shape, dtype=numpy.float32)
        asmos_revised[:] = _TARGET_NODATA
        asmos_revised[valid_mask] = asmos_interm[valid_mask]
        asmos_revised[exceeded_mask] = afl[exceeded_mask]

        notexceeded_mask = (valid_mask & (asmos_interm <= afl))
        amov[notexceeded_mask] = 0.

        if return_type == 'asmos_revised':
            return asmos_revised
        elif return_type == 'amov':
            return amov
    return _distribute_water


def calc_available_water_for_transpiration(asmos, awilt, adep):
    """Calculate water available for transpiration in one soil layer.

    The water available for transpiration is the amount of water in the soil
    layer minus the wilting point of the soil layer.

    Parameters:
        asmos (numpy.ndarray): derived, interim moisture in the soil layer
        awilt (numpy.ndarray): derived, wilting point of the soil layer
        adep (numpy.ndarray): parameter, depth of the soil layer in cm

    Returns:
        avw, available water for transpiration
    """
    valid_mask = (
        (asmos != _TARGET_NODATA) &
        (awilt != _TARGET_NODATA) &
        (adep != _IC_NODATA))
    avw = numpy.empty(asmos.shape, dtype=numpy.float32)
    avw[:] = _TARGET_NODATA
    avw[valid_mask] = numpy.maximum(
        asmos[valid_mask] - awilt[valid_mask] * adep[valid_mask], 0.)
    return avw


def revise_potential_transpiration(trap, tot):
    """Revise potential transpiration according to water available.

    Total potential transpiration, trap, is revised to be less than or equal
    to total water available for transpiration, tot. Total water available
    for transpiration is the sum of available water per soil layer.
    Line 241, H2olos.f

    Parameters:
        trap (numpy.ndarray): derived, potential transpiration water losses
        tot (numpy.ndarray): derived, total soil water available for
            transpiration

    Returns:
        trap_revised, revised total potential transpiration
    """
    valid_mask = (
        (trap != _TARGET_NODATA) &
        (tot != _TARGET_NODATA))
    trap_revised = numpy.empty(trap.shape, dtype=numpy.float32)
    trap_revised[:] = _TARGET_NODATA
    trap_revised[valid_mask] = numpy.minimum(trap[valid_mask], tot[valid_mask])
    return trap_revised


def remove_transpiration(return_type):
    """Remove water from a soil layer via transpiration by plants.

    Transpiration from one soil layer is apportioned from total potential
    transpiration, trap, according to the available water for transpiration in
    this soil layer. Lines 218-294, H2olos.f

    Parameters:
        return_type (string): flag indicating whether avinj (water in this soil
            layer available to plants for growth) or asmos (total water in this
            soil layer) should be returned

    Returns:
        the function `_remove_transpiration`
    """
    def _remove_transpiration(asmos, awilt, adep, trap, awwt, tot2):
        """Remove water from a soil layer via transpiration by plants.

        Parameters:
            asmos (numpy.ndarray): derived, interim moisture in this soil layer
                after additions from current month precipitation
            awilt (numpy.ndarray): derived, wilting point of this soil layer
            adep (numpy.ndarray): parameter, depth of this soil layer in cm
            trap (numpy.ndarray): derived, total potential transpiration
                across all soil layers accessible by plant roots
            awwt (numpy.ndarray): derived, water available for transpiration
                in this soil layer weighted by transpiration depth distribution
                parameter
            tot2 (numpy.ndarray): derived, the sum of weighted water available
                for transpiration across soil layers

        Returns:
            avinj, water available to plants for growth in this layer after
                losses to transpiration, if return type is 'avinj'
            asmos_revised, total water in this layer after losses to
                transpiration, if return type is 'asmos'
        """
        valid_mask = (
            (asmos != _TARGET_NODATA) &
            (awilt != _TARGET_NODATA) &
            (adep != _IC_NODATA) &
            (trap != _TARGET_NODATA) &
            (awwt != _TARGET_NODATA) &
            (tot2 != _TARGET_NODATA))

        avinj = numpy.empty(asmos.shape, dtype=numpy.float32)
        avinj[:] = _TARGET_NODATA
        avinj[valid_mask] = numpy.maximum(
            asmos[valid_mask] - awilt[valid_mask] * adep[valid_mask], 0.)

        transpire_mask = (valid_mask & (tot2 > 0))
        transpiration_loss = numpy.zeros(asmos.shape, dtype=numpy.float32)
        transpiration_loss[transpire_mask] = numpy.minimum(
            (trap[transpire_mask] *
                awwt[transpire_mask]) / tot2[transpire_mask],
            avinj[transpire_mask])
        avinj[valid_mask] = avinj[valid_mask] - transpiration_loss[valid_mask]

        asmos_revised = numpy.empty(asmos.shape, dtype=numpy.float32)
        asmos_revised[:] = _TARGET_NODATA
        asmos_revised[valid_mask] = (
            asmos[valid_mask] - transpiration_loss[valid_mask])

        if return_type == 'avinj':
            return avinj
        elif return_type == 'asmos':
            return asmos_revised
    return _remove_transpiration


def calc_relative_water_content_lyr_1(asmos_1, adep_1, awilt_1, afiel_1):
    """Calculate the relative water content of soil layer 1.

    The relative water content of soil layer 1, prior to any evaporation losses
    from soil layer 1, is used to estimate water available for evaporation
    from soil layer 1. Line 280, H2olos.f

    Parameters:
        asmos_1 (numpy.ndarray): derived, interim moisture in soil layer 1
            after losses to transpiration
        adep_1 (numpy.ndarray): parameter, depth of soil layer 1 in cm
        awilt_1 (numpy.ndarray): derived, wilting point of soil layer 1
        afiel_1 (numpy.ndarray): derived, field capacity of soil layer 1

    Returns:
        rwcf_1, relative water content of soil layer 1
    """
    valid_mask = (
        (asmos_1 != _TARGET_NODATA) &
        (adep_1 != _IC_NODATA) &
        (awilt_1 != _TARGET_NODATA) &
        (afiel_1 != _TARGET_NODATA))
    rwcf_1 = numpy.empty(asmos_1.shape, dtype=numpy.float32)
    rwcf_1[valid_mask] = (
        (asmos_1[valid_mask] / adep_1[valid_mask] - awilt_1[valid_mask]) /
        (afiel_1[valid_mask] - awilt_1[valid_mask]))
    return rwcf_1


def calc_evaporation_loss(rwcf_1, pevp, absevap, asmos_1, awilt_1, adep_1):
    """Calculate evaporation from soil layer 1.

    Some moisture is lost from soil layer 1 (i.e., the top soil layer) to
    evaporation, separate from surface evaporation and transpiration by plants.
    This amount is calculated from potential soil evaporation, which was
    calculated from potential evapotranspiration prior to allocation of water
    to soil layers. It is restricted to be less than or equal to water
    available in this soil layer.

    Parameters:
        rwcf_1 (numpy.ndarray): derived, relative water content of soil layer 1
        pevp (numpy.ndarray): derived, potential evaporation from soil layer 1
        absevap (numpy.ndarray): derived, bare soil evaporation
        asmos_1 (numpy.ndarray): derived, interim moisture in soil layer 1
        awilt_1 (numpy.ndarray): derived, wilting point of soil layer 1
        adep_1 (numpy.ndarray): parameter, depth of soil layer 1 in cm

    Returns:
        evlos, moisture evaporated from soil layer 1
    """
    valid_mask = (
        (rwcf_1 != _TARGET_NODATA) &
        (pevp != _TARGET_NODATA) &
        (absevap != _TARGET_NODATA) &
        (asmos_1 != _TARGET_NODATA) &
        (awilt_1 != _TARGET_NODATA) &
        (adep_1 != _IC_NODATA))
    evmt = numpy.empty(rwcf_1.shape, dtype=numpy.float32)
    evmt[:] = _TARGET_NODATA
    evmt[valid_mask] = numpy.maximum(
        (rwcf_1[valid_mask] - 0.25) / (1 - 0.25), 0.01)
    evlos = numpy.empty(rwcf_1.shape, dtype=numpy.float32)
    evlos[:] = _TARGET_NODATA
    evlos[valid_mask] = numpy.minimum(
        evmt[valid_mask] * pevp[valid_mask] * absevap[valid_mask] * 0.1,
        numpy.maximum(
            asmos_1[valid_mask] - awilt_1[valid_mask] *
            adep_1[valid_mask], 0.))
    return evlos


def _soil_water(
        aligned_inputs, site_param_table, veg_trait_table, current_month,
        month_index, prev_sv_reg, sv_reg, pp_reg, month_reg, pft_id_set):
    """Allocate precipitation to runoff, transpiration, and soil moisture.

    Simulate snowfall and account for evaporation and melting of the snow pack.
    Allocate the flow of precipitation through interception by plants,
    runoff and infiltration into the soil, percolation through the soil, and
    transpiration by plants. Update soil moisture in each soil layer.
    Estimate avh2o_1 for each PFT (water available to the PFT for growth),
    avh2o_3 (water in first two soil layers), and amov_2 (indicating whether
    movement of water occurs from the second soil layer, used in
    decomposition).

    Parameters:
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including precipitation, temperature,
            plant functional type composition, and site spatial index
        site_param_table (dict): map of site spatial indices to dictionaries
            containing site parameters
        veg_trait_table (dict): map of pft id to dictionaries containing
            plant functional type parameters, including nlaypg, number of soil
            layers access by plant roots
        current_month (int): month of the year, such that current_month=1
            indicates January
        month_index (int): month of the simulation, such that month_index=1
            indicates month 1 of the simulation
        prev_sv_reg (dict): map of key, path pairs giving paths to state
            variables for the previous month
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month
        pp_reg (dict): map of key, path pairs giving persistent parameters
            including field capacity of each soil layer
        month_reg (dict): map of key, path pairs giving paths to intermediate
            calculated values that are shared between submodels
        pft_id_set (set): set of integers identifying plant functional types

    Modifies:
        the raster indicated by `sv_reg['snow_path']`, current snowpack
        the raster indicated by `sv_reg['snlq_path']`, current liquid in snow
        the raster indicated by `sv_reg['asmos_<lyr>_path']`, soil moisture
            content, for each soil layer accessible by roots of any plant
            functional type
        the raster indicated by `month_reg['amov_2']`, movement of soil
            moisture out of soil layer 2
        the raster indicated by `sv_reg['avh2o_1_<PFT>_path']`, soil moisture
            available for growth, for each plant functional type (PFT)
        the raster indicated by `sv_reg['avh2o_3_path']`, available water in
            the top two soil layers

    Returns:
        None
    """
    def calc_avg_temp(max_temp, min_temp):
        """Calculate average temperature from maximum and minimum temp."""
        valid_mask = (
            (~numpy.isclose(max_temp, max_temp_nodata)) &
            (~numpy.isclose(min_temp, min_temp_nodata)))
        tave = numpy.empty(max_temp.shape, dtype=numpy.float32)
        tave[:] = _IC_NODATA
        tave[valid_mask] = (max_temp[valid_mask] + min_temp[valid_mask]) / 2.
        return tave

    def calc_surface_litter_biomass(strucc_1, metabc_1):
        """Calculate biomass in surface litter."""
        valid_mask = (
            (~numpy.isclose(strucc_1, _SV_NODATA)) &
            (~numpy.isclose(metabc_1, _SV_NODATA)))
        alit = numpy.empty(strucc_1.shape, dtype=numpy.float32)
        alit[:] = _TARGET_NODATA
        alit[valid_mask] = (strucc_1[valid_mask] + metabc_1[valid_mask]) * 2.5
        alit = numpy.minimum(alit, 400)
        return alit

    max_temp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['max_temp_{}'.format(current_month)])['nodata'][0]
    min_temp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['min_temp_{}'.format(current_month)])['nodata'][0]

    # get max number of soil layers to track
    nlaypg_max = int(max(
        veg_trait_table[i]['nlaypg'] for i in veg_trait_table.iterkeys()))

    # temporary intermediate rasters for soil water submodel
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for val in [
            'tave', 'current_moisture_inputs', 'modified_moisture_inputs',
            'pet_rem', 'alit', 'sum_aglivc', 'sum_stdedc', 'sum_tgprod',
            'aliv', 'sd', 'absevap', 'evap_losses', 'trap', 'trap_revised',
            'pevp', 'tot', 'tot2', 'rwcf_1', 'evlos', 'avinj_interim_1']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))
    for val in ['asmos_interim', 'avw', 'awwt', 'avinj']:
        for lyr in xrange(1, nlaypg_max + 1):
            val_lyr = '{}_{}'.format(val, lyr)
            temp_val_dict[val_lyr] = os.path.join(
                temp_dir, '{}.tif'.format(val_lyr))
    # PFT-level temporary calculated values
    for pft_i in pft_id_set:
        for val in ['tgprod_weighted', 'sum_avinj']:
            temp_val_dict['{}_{}'.format(val, pft_i)] = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))

    param_val_dict = {}
    for val in ['fracro', 'precro', 'fwloss_1', 'fwloss_2']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (aligned_inputs['site_index'], 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)
    for val in ['adep', 'awtl']:
        for lyr in xrange(1, nlaypg_max + 1):
            val_lyr = '{}_{}'.format(val, lyr)
            target_path = os.path.join(temp_dir, '{}.tif'.format(val_lyr))
            param_val_dict[val_lyr] = target_path
            site_to_val = dict(
                [(site_code, float(table[val_lyr])) for
                    (site_code, table) in site_param_table.iteritems()])
            pygeoprocessing.reclassify_raster(
                (aligned_inputs['site_index'], 1), site_to_val, target_path,
                gdal.GDT_Float32, _IC_NODATA)

    # calculate canopy and litter cover that influence moisture inputs
    # calculate biomass in surface litter
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            prev_sv_reg['strucc_1_path'], prev_sv_reg['metabc_1_path']]],
        calc_surface_litter_biomass, temp_val_dict['alit'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate the sum of aglivc (standing live biomass) and stdedc
    # (standing dead biomass) across PFTs, weighted by % cover of each PFT
    for sv in ['aglivc', 'stdedc']:
        weighted_sum_path = temp_val_dict['sum_{}'.format(sv)]
        weighted_state_variable_sum(
            sv, prev_sv_reg, aligned_inputs, pft_id_set, weighted_sum_path)

    # calculate the weighted sum of tgprod, potential production, across PFTs
    weighted_path_list = []
    for pft_i in pft_id_set:
        target_path = temp_val_dict['tgprod_weighted_{}'.format(pft_i)]
        if month_reg['tgprod_{}'.format(pft_i)]:
            pft_nodata = pygeoprocessing.get_raster_info(
                aligned_inputs['pft_{}'.format(pft_i)])['nodata'][0]
            raster_multiplication(
                month_reg['tgprod_{}'.format(pft_i)], _TARGET_NODATA,
                aligned_inputs['pft_{}'.format(pft_i)], pft_nodata,
                target_path, _TARGET_NODATA)
            weighted_path_list.append(target_path)
    if weighted_path_list:
        raster_list_sum(
            weighted_path_list, _TARGET_NODATA,
            temp_val_dict['sum_tgprod'], _TARGET_NODATA, nodata_remove=True)
    else:  # no potential production occurs this month, so tgprod = 0
        pygeoprocessing.new_raster_from_base(
            temp_val_dict['sum_aglivc'], temp_val_dict['sum_tgprod'],
            gdal.GDT_Float32, [_TARGET_NODATA], fill_value_list=[0.])

    # calculate average temperature
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            aligned_inputs['max_temp_{}'.format(current_month)],
            aligned_inputs['min_temp_{}'.format(current_month)]]],
        calc_avg_temp, temp_val_dict['tave'], gdal.GDT_Float32, _IC_NODATA)

    # calculate aboveground live biomass
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['sum_aglivc'], temp_val_dict['sum_tgprod']]],
        _calc_aboveground_live_biomass, temp_val_dict['aliv'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate total standing biomass
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['aliv'], temp_val_dict['sum_stdedc']]],
        _calc_standing_biomass, temp_val_dict['sd'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # modify standing snow, liquid in snow, return moisture inputs after snow
    _snow(
        aligned_inputs['site_index'], site_param_table,
        aligned_inputs['precip_{}'.format(month_index)],
        temp_val_dict['tave'],
        aligned_inputs['max_temp_{}'.format(current_month)],
        aligned_inputs['min_temp_{}'.format(current_month)],
        prev_sv_reg['snow_path'], prev_sv_reg['snlq_path'],
        current_month, month_reg['snowmelt'], sv_reg['snow_path'],
        sv_reg['snlq_path'], temp_val_dict['modified_moisture_inputs'],
        temp_val_dict['pet_rem'])

    # remove runoff and surface evaporation from moisture inputs
    shutil.copyfile(
        temp_val_dict['modified_moisture_inputs'],
        temp_val_dict['current_moisture_inputs'])
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['current_moisture_inputs'],
            param_val_dict['fracro'], param_val_dict['precro'],
            sv_reg['snow_path'], temp_val_dict['alit'],
            temp_val_dict['sd'], param_val_dict['fwloss_1'],
            param_val_dict['fwloss_2'], temp_val_dict['pet_rem']]],
        subtract_surface_losses('inputs_after_surface'),
        temp_val_dict['modified_moisture_inputs'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate bare soil evaporation
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['current_moisture_inputs'],
            param_val_dict['fracro'], param_val_dict['precro'],
            sv_reg['snow_path'], temp_val_dict['alit'],
            temp_val_dict['sd'], param_val_dict['fwloss_1'],
            param_val_dict['fwloss_2'], temp_val_dict['pet_rem']]],
        subtract_surface_losses('absevap'),
        temp_val_dict['absevap'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate total losses to surface evaporation
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['current_moisture_inputs'],
            param_val_dict['fracro'], param_val_dict['precro'],
            sv_reg['snow_path'], temp_val_dict['alit'],
            temp_val_dict['sd'], param_val_dict['fwloss_1'],
            param_val_dict['fwloss_2'], temp_val_dict['pet_rem']]],
        subtract_surface_losses('evap_losses'),
        temp_val_dict['evap_losses'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # remove losses due to initial transpiration from water inputs
    shutil.copyfile(
        temp_val_dict['modified_moisture_inputs'],
        temp_val_dict['current_moisture_inputs'])
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['pet_rem'], temp_val_dict['evap_losses'],
            temp_val_dict['tave'], temp_val_dict['aliv'],
            temp_val_dict['current_moisture_inputs']]],
        calc_potential_transpiration('modified_moisture_inputs'),
        temp_val_dict['modified_moisture_inputs'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate potential transpiration
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['pet_rem'], temp_val_dict['evap_losses'],
            temp_val_dict['tave'], temp_val_dict['aliv'],
            temp_val_dict['current_moisture_inputs']]],
        calc_potential_transpiration('trap'), temp_val_dict['trap'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # calculate potential evaporation from top soil layer
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['pet_rem'], temp_val_dict['evap_losses'],
            temp_val_dict['tave'], temp_val_dict['aliv'],
            temp_val_dict['current_moisture_inputs']]],
        calc_potential_transpiration('pevp'), temp_val_dict['pevp'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # distribute water to each layer
    for lyr in xrange(1, nlaypg_max + 1):
        shutil.copyfile(
            temp_val_dict['modified_moisture_inputs'],
            temp_val_dict['current_moisture_inputs'])
        # revise moisture content of this soil layer
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                param_val_dict['adep_{}'.format(lyr)],
                pp_reg['afiel_{}_path'.format(lyr)],
                prev_sv_reg['asmos_{}_path'.format(lyr)],
                temp_val_dict['current_moisture_inputs']]],
            distribute_water_to_soil_layer('asmos_revised'),
            temp_val_dict['asmos_interim_{}'.format(lyr)],
            gdal.GDT_Float32, _TARGET_NODATA)
        # calculate soil moisture moving to next layer
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                param_val_dict['adep_{}'.format(lyr)],
                pp_reg['afiel_{}_path'.format(lyr)],
                prev_sv_reg['asmos_{}_path'.format(lyr)],
                temp_val_dict['current_moisture_inputs']]],
            distribute_water_to_soil_layer('amov'),
            temp_val_dict['modified_moisture_inputs'],
            gdal.GDT_Float32, _TARGET_NODATA)
        if lyr == 2:
            # moisture leaving soil layer 2 is needed in decomposition submodel
            shutil.copyfile(
                temp_val_dict['modified_moisture_inputs'], month_reg['amov_2'])

    # calculate available water for transpiration
    avw_list = []
    for lyr in xrange(1, nlaypg_max + 1):
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['asmos_interim_{}'.format(lyr)],
                pp_reg['awilt_{}_path'.format(lyr)],
                param_val_dict['adep_{}'.format(lyr)]]],
            calc_available_water_for_transpiration,
            temp_val_dict['avw_{}'.format(lyr)], gdal.GDT_Float32,
            _TARGET_NODATA)
        avw_list.append(temp_val_dict['avw_{}'.format(lyr)])
    # total water available for transpiration
    raster_list_sum(
        avw_list, _TARGET_NODATA, temp_val_dict['tot'], _TARGET_NODATA)

    # calculate water available for transpiration weighted by transpiration
    # depth for that soil layer
    awwt_list = []
    for lyr in xrange(1, nlaypg_max + 1):
        raster_multiplication(
            temp_val_dict['avw_{}'.format(lyr)], _TARGET_NODATA,
            param_val_dict['awtl_{}'.format(lyr)], _IC_NODATA,
            temp_val_dict['awwt_{}'.format(lyr)], _TARGET_NODATA)
        awwt_list.append(temp_val_dict['awwt_{}'.format(lyr)])
    # total weighted available water for transpiration
    raster_list_sum(
        awwt_list, _TARGET_NODATA, temp_val_dict['tot2'], _TARGET_NODATA)

    # revise total potential transpiration
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [temp_val_dict['trap'], temp_val_dict['tot']]],
        revise_potential_transpiration, temp_val_dict['trap_revised'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # remove water via transpiration
    for lyr in xrange(1, nlaypg_max + 1):
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['asmos_interim_{}'.format(lyr)],
                pp_reg['awilt_{}_path'.format(lyr)],
                param_val_dict['adep_{}'.format(lyr)],
                temp_val_dict['trap_revised'],
                temp_val_dict['awwt_{}'.format(lyr)], temp_val_dict['tot2']]],
            remove_transpiration('avinj'),
            temp_val_dict['avinj_{}'.format(lyr)], gdal.GDT_Float32,
            _TARGET_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['asmos_interim_{}'.format(lyr)],
                pp_reg['awilt_{}_path'.format(lyr)],
                param_val_dict['adep_{}'.format(lyr)],
                temp_val_dict['trap_revised'],
                temp_val_dict['awwt_{}'.format(lyr)], temp_val_dict['tot2']]],
            remove_transpiration('asmos'), sv_reg['asmos_{}_path'.format(lyr)],
            gdal.GDT_Float32, _TARGET_NODATA)

    # relative water content of soil layer 1
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            sv_reg['asmos_1_path'], param_val_dict['adep_1'],
            pp_reg['awilt_1_path'], pp_reg['afiel_1_path']]],
        calc_relative_water_content_lyr_1, temp_val_dict['rwcf_1'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # evaporation from soil layer 1
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['rwcf_1'], temp_val_dict['pevp'],
            temp_val_dict['absevap'], sv_reg['asmos_1_path'],
            pp_reg['awilt_1_path'], param_val_dict['adep_1']]],
        calc_evaporation_loss, temp_val_dict['evlos'],
        gdal.GDT_Float32, _TARGET_NODATA)

    # remove evaporation from total moisture in soil layer 1
    shutil.copyfile(sv_reg['asmos_1_path'], temp_val_dict['asmos_interim_1'])
    raster_difference(
        temp_val_dict['asmos_interim_1'], _TARGET_NODATA,
        temp_val_dict['evlos'], _TARGET_NODATA, sv_reg['asmos_1_path'],
        _TARGET_NODATA)

    # remove evaporation from moisture available to plants in soil layer 1
    shutil.copyfile(temp_val_dict['avinj_1'], temp_val_dict['avinj_interim_1'])
    raster_difference(
        temp_val_dict['avinj_interim_1'], _TARGET_NODATA,
        temp_val_dict['evlos'], _TARGET_NODATA, temp_val_dict['avinj_1'],
        _TARGET_NODATA)

    # calculate avh2o_1, soil water available for growth, for each PFT
    for pft_i in pft_id_set:
        pft_nodata = pygeoprocessing.get_raster_info(
            aligned_inputs['pft_{}'.format(pft_i)])['nodata'][0]
        soil_layers_accessible = [
            temp_val_dict['avinj_{}'.format(lyr)] for lyr in
            xrange(1, int(veg_trait_table[pft_i]['nlaypg']) + 1)]
        raster_list_sum(
            soil_layers_accessible, _TARGET_NODATA,
            temp_val_dict['sum_avinj_{}'.format(pft_i)],
            _TARGET_NODATA, nodata_remove=True)
        raster_multiplication(
            temp_val_dict['sum_avinj_{}'.format(pft_i)], _TARGET_NODATA,
            aligned_inputs['pft_{}'.format(pft_i)], pft_nodata,
            sv_reg['avh2o_1_{}_path'.format(pft_i)], _SV_NODATA)

    # calculate avh2o_3, moisture in top two soil layers
    soil_layers_to_sum = [
        temp_val_dict['avinj_{}'.format(lyr)] for lyr in [1, 2]]
    raster_list_sum(
        soil_layers_to_sum, _TARGET_NODATA, sv_reg['avh2o_3_path'],
        _SV_NODATA, nodata_remove=False)

    # set correct nodata value for all revised asmos rasters
    for lyr in xrange(1, nlaypg_max + 1):
        reclassify_nodata(sv_reg['asmos_{}_path'.format(lyr)], _SV_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def _monthly_N_fixation(
        aligned_inputs, month_index, site_param_table, year_reg, prev_sv_reg,
        sv_reg):
    """Add monthly atmospheric nitrogen fixation to surface mineral N.

    Atmospheric N fixation for the month is calculated from annual N
    deposition, calculated once per year, according to the ratio of monthly
    precipitation to annual precipitation.  Total N fixed in this month is
    added to the surface mineral N pool.  Lines 193-205, Simsom.f

    Parameters:
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including precipitation and site index
            path
        month_index (int): month of the simulation, such that month_index=13
            indicates month 13 of the simulation
        site_param_table (dict): map of site spatial indices to dictionaries
            containing site parameters
        year_reg (dict): map of key, path pairs giving paths to rasters that
            are modified once per year, including annual precipitation and base
            N deposition
        prev_sv_reg (dict): map of key, path pairs giving paths to state
            variables for the previous month
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month, including minerl_1_1, mineral N in the
            surface layer

    Modifies:
        the raster indicated by sv_reg['minerl_1_1_path']

    Returns:
        none
    """
    def calc_N_fixation(
            precip, annual_precip, baseNdep, epnfs_2, prev_minerl_1_1):
        """Add monthly N fixation to surface mineral N pool.

        Monthly N fixation is calculated from annual N deposition according to
        the ratio of monthly precipitation to annual precipitation.

        Parameters:
            precip (numpy.ndarray): input, monthly precipitation
            annual_precip (numpy.ndarray): derived, annual precipitation
            baseNdep (numpy.ndarray): derived, annual atmospheric N deposition
            epnfs_2 (numpy.ndarray): parameter, intercept of regression
                predicting N deposition from annual precipitation
            prev_minerl_1_1 (numpy.ndarray): state variable, mineral N in the
                surface layer in previous month

        Returns:
            minerl_1_1, updated mineral N in the surface layer
        """
        valid_mask = (
            (~numpy.isclose(precip, precip_nodata)) &
            (annual_precip != _TARGET_NODATA) &
            (baseNdep != _TARGET_NODATA) &
            (epnfs_2 != _IC_NODATA) &
            (~numpy.isclose(prev_minerl_1_1, _SV_NODATA)))
        wdfxm = numpy.zeros(precip.shape, dtype=numpy.float32)
        wdfxm[valid_mask] = (
            baseNdep[valid_mask] *
            (precip[valid_mask] / annual_precip[valid_mask]) +
            epnfs_2[valid_mask] *
            numpy.minimum(annual_precip[valid_mask], 100.) *
            (precip[valid_mask] / annual_precip[valid_mask]))

        minerl_1_1 = numpy.empty(precip.shape, dtype=numpy.float32)
        minerl_1_1[:] = _SV_NODATA
        minerl_1_1[valid_mask] = (
            prev_minerl_1_1[valid_mask] + wdfxm[valid_mask])
        return minerl_1_1

    precip_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['precip_{}'.format(month_index)])['nodata'][0]

    # temporary intermediate rasters for calculating available nutrient
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    target_path = os.path.join(temp_dir, 'epnfs_2.tif')
    param_val_dict = {'epnfs_2': target_path}
    site_to_val = dict(
        [(site_code, float(table['epnfs_2'])) for
            (site_code, table) in site_param_table.iteritems()])
    pygeoprocessing.reclassify_raster(
        (aligned_inputs['site_index'], 1), site_to_val, target_path,
        gdal.GDT_Float32, _IC_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            aligned_inputs['precip_{}'.format(month_index)],
            year_reg['annual_precip_path'], year_reg['baseNdep_path'],
            param_val_dict['epnfs_2'], prev_sv_reg['minerl_1_1_path']]],
        calc_N_fixation, sv_reg['minerl_1_1_path'], gdal.GDT_Float32,
        _SV_NODATA)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def calc_anerb(rprpet, pevap, drain, aneref_1, aneref_2, aneref_3):
    """Calculate the effect of soil anaerobic conditions on decomposition.

    The impact of soil anaerobic conditions on decomposition is calculated from
    soil moisture and reference evapotranspiration. Anerob.f.

    Parameters:
        rprpet (numpy.ndarray): derived, ratio of precipitation or snowmelt to
            reference evapotranspiration
        pevap (numpy.ndarray): derived, reference evapotranspiration
        drain (numpy.ndarray): parameter, the fraction of excess water lost by
            drainage. Indicates whether a soil is sensitive for anaerobiosis
            (drain = 0) or not (drain = 1)
        aneref_1 (numpy.ndarray): parameter, value of rprpet below which there
            is no negative impact of soil anaerobic conditions on decomposition
        aneref_2 (numpy.ndarray): parameter, value of rprpet above which there
            is maximum negative impact of soil anaerobic conditions on
            decomposition
        aneref_3 (numpy.ndarray): parameter, minimum value of the impact of
            soil anaerobic conditions on decomposition

    Returns:
        anerb, the effect of soil anaerobic conditions on decomposition
    """
    valid_mask = (
        (rprpet != _TARGET_NODATA) &
        (pevap != _TARGET_NODATA) &
        (drain != _IC_NODATA) &
        (aneref_1 != _IC_NODATA) &
        (aneref_2 != _IC_NODATA) &
        (aneref_3 != _IC_NODATA))

    xh2o = numpy.empty(rprpet.shape, dtype=numpy.float32)
    xh2o[:] = _TARGET_NODATA
    xh2o[valid_mask] = (
        (rprpet[valid_mask] - aneref_1[valid_mask]) * pevap[valid_mask] *
        (1. - drain[valid_mask]))

    anerb = numpy.empty(rprpet.shape, dtype=numpy.float32)
    anerb[:] = _TARGET_NODATA
    anerb[valid_mask] = 1.

    high_rprpet_mask = (valid_mask & (rprpet > aneref_1) & (xh2o > 0))
    anerb[high_rprpet_mask] = numpy.maximum(
        1. + (1. - aneref_3[high_rprpet_mask]) /
        (aneref_1[high_rprpet_mask] - aneref_2[high_rprpet_mask]) *
        (aneref_1[high_rprpet_mask] +
            (xh2o[high_rprpet_mask] / pevap[high_rprpet_mask]) -
            aneref_1[high_rprpet_mask]),
        aneref_3[high_rprpet_mask])
    return anerb


def esched(return_type):
    """Calculate flow of an element accompanying decomposition of C.

    Calculate the movement of one element (N or P) as C decomposes from one
    state variable (the donating stock, or box A) to another state variable
    (the receiving stock, or box B).  Esched.f

    Parameters:
        return_type (string): flag indicating whether to return material
            leaving box A, material arriving in box B, or material flowing
            into or out of the mineral pool

    Returns:
        the function `_esched`
    """
    def _esched(cflow, tca, rcetob, anps, labile):
        """Calculate the flow of one element (iel) to accompany decomp of C.

        This is a transcription of Esched.f: "Schedule N, P, or S flow and
        associated mineralization or immobilization flow for decomposition
        from Box A to Box B."
        If there is enough of iel (N or P) in the donating stock to satisfy
        the required ratio, that material flows from the donating stock to
        the receiving stock and whatever iel is leftover goes to mineral
        pool. If there is not enough iel to satisfy the required ratio, iel
        is drawn from the mineral pool to satisfy the ratio; if there is
        not enough iel in the mineral pool, the material does not leave the
        donating stock.

        Parameters:
            cflow (numpy.ndarray): derived, total C that is decomposing from
                box A to box B
            tca (numpy.ndarray): state variable, C in donating stock, i.e.
                box A
            rcetob (numpy.ndarray): derived, required ratio of C/iel in the
                receiving stock
            anps (numpy.ndarray): state variable, iel (N or P) in the donating
                stock
            labile (numpy.ndarray): state variable, mineral iel (N or P)

        Returns:
            material_leaving_a, the amount of material leaving box A, if
                return_type is 'material_leaving_a'
            material_arriving_b, the amount of material arriving in box B,
                if return_type is 'material_arriving_b'
            mnrflo, flow to or from mineral pool, if return_type is
                'mineral_flow'
        """
        valid_mask = (
            (cflow != _IC_NODATA) &
            (~numpy.isclose(tca, _SV_NODATA)) &
            (rcetob != _TARGET_NODATA) &
            (~numpy.isclose(anps, _SV_NODATA)) &
            (~numpy.isclose(labile, _SV_NODATA)))
        outofa = numpy.empty(cflow.shape, dtype=numpy.float32)
        outofa[:] = _IC_NODATA
        outofa[valid_mask] = (
            anps[valid_mask] * (cflow[valid_mask] / tca[valid_mask]))

        immobil_ratio = numpy.zeros(cflow.shape)
        immobil_ratio[valid_mask] = (
            cflow[valid_mask] / outofa[valid_mask])

        immflo = numpy.zeros(cflow.shape)
        immflo[valid_mask] = (
            cflow[valid_mask] / rcetob[valid_mask] - outofa[valid_mask])

        labile_supply = numpy.zeros(cflow.shape)
        labile_supply[valid_mask] = labile[valid_mask] - immflo[valid_mask]

        atob = numpy.zeros(cflow.shape)
        atob[valid_mask] = cflow[valid_mask] / rcetob[valid_mask]

        # immobilization
        immobilization_mask = (
            (immobil_ratio > rcetob) &
            (labile_supply > 0) &
            valid_mask)
        # mineralization
        mineralization_mask = (
            (immobil_ratio <= rcetob) &
            valid_mask)
        # no movement
        no_movt_mask = (
            (immobil_ratio > rcetob) &
            (labile_supply <= 0) &
            valid_mask)

        material_leaving_a = numpy.empty(cflow.shape, dtype=numpy.float32)
        material_leaving_a[:] = _IC_NODATA
        material_arriving_b = numpy.empty(cflow.shape, dtype=numpy.float32)
        material_arriving_b[:] = _IC_NODATA
        mnrflo = numpy.empty(cflow.shape, dtype=numpy.float32)
        mnrflo[:] = _IC_NODATA

        material_leaving_a[immobilization_mask] = (
            outofa[immobilization_mask])
        material_arriving_b[immobilization_mask] = (
            outofa[immobilization_mask] + immflo[immobilization_mask])
        mnrflo[immobilization_mask] = -immflo[immobilization_mask]

        material_leaving_a[mineralization_mask] = outofa[mineralization_mask]
        material_arriving_b[mineralization_mask] = atob[mineralization_mask]
        mnrflo[mineralization_mask] = (
            outofa[mineralization_mask] - atob[mineralization_mask])

        material_leaving_a[no_movt_mask] = 0.
        material_arriving_b[no_movt_mask] = 0.
        mnrflo[no_movt_mask] = 0.

        if return_type == 'material_leaving_a':
            return material_leaving_a
        elif return_type == 'material_arriving_b':
            return material_arriving_b
        elif return_type == 'mineral_flow':
            return mnrflo
    return _esched


def fsfunc(minerl_1_2, sorpmx, pslsrb):
    """Calculate the fraction of mineral P that is in solution.

    The fraction of P in solution is influenced by two soil properties:
    the maximum sorption potential of the soil and sorption affinity.

    Parameters:
        minerl_1_2 (numpy.ndarray): state variable, mineral P in top layer
        sorpmx (numpy.ndarray): parameter, maximum P sorption potential
        pslsrb (numpy.ndarray): parameter, slope term which controls the
            fraction of mineral P that is labile

    Returns:
        fsol, fraction of P in solution
    """
    valid_mask = (
        (~numpy.isclose(minerl_1_2, _SV_NODATA)) &
        (sorpmx != _IC_NODATA) &
        (pslsrb != _IC_NODATA))
    c_ar = numpy.zeros(minerl_1_2.shape, dtype=numpy.float32)
    c_ar[valid_mask] = (
        sorpmx[valid_mask] * (2.0 - pslsrb[valid_mask]) / 2.)
    b_ar = numpy.zeros(minerl_1_2.shape, dtype=numpy.float32)
    b_ar[valid_mask] = (
        sorpmx[valid_mask] - minerl_1_2[valid_mask] + c_ar[valid_mask])

    sq_ar = numpy.zeros(minerl_1_2.shape, dtype=numpy.float32)
    sq_ar[valid_mask] = (
        b_ar[valid_mask] * b_ar[valid_mask] + 4. * c_ar[valid_mask] *
        minerl_1_2[valid_mask])
    sqrt_ar = numpy.zeros(minerl_1_2.shape, dtype=numpy.float32)
    sqrt_ar[valid_mask] = numpy.sqrt(sq_ar[valid_mask])

    labile = numpy.zeros(minerl_1_2.shape, dtype=numpy.float32)
    labile[valid_mask] = (-b_ar[valid_mask] + sqrt_ar[valid_mask]) / 2.

    fsol = numpy.empty(minerl_1_2.shape, dtype=numpy.float32)
    fsol[:] = _TARGET_NODATA
    fsol[valid_mask] = labile[valid_mask] / minerl_1_2[valid_mask]
    return fsol


def calc_surface_som2_ratio(
        som1c_1, som1e_1_iel, rad1p_1_iel, rad1p_2_iel, rad1p_3_iel,
        pcemic1_2_iel):
    """Calculate the required C/iel ratio for material entering surface SOM2.

    The C/iel ratio of material decomposing from surface SOM1 into surface SOM2
    fluctuates with each decomposition time step according to the current C/iel
    content of SOM1.

    Parameters:
        som1c_1 (numpy.ndarray): state variable, C in surface SOM1
        som1e_1_iel (numpy.ndarray): state variable, iel in surface SOM1
        rad1p_1_iel (numpy.ndarray): parameter, intercept term
        rad1p_2_iel (numpy.ndarray): parameter, slope term
        rad1p_3_iel (numpy.ndarray): parameter, minimum allowable C/iel for
            addition term
        pcemic1_2_iel (numpy.ndarray): parameter, minimum C/iel ratio

    Returns:
        rceto2_surface, required C/iel ratio of material entering surface SOM2
    """
    valid_mask = (
        (~numpy.isclose(som1c_1, _SV_NODATA)) &
        (~numpy.isclose(som1e_1_iel, _SV_NODATA)) &
        (rad1p_1_iel != _IC_NODATA) &
        (rad1p_2_iel != _IC_NODATA) &
        (pcemic1_2_iel != _IC_NODATA) &
        (rad1p_3_iel != _IC_NODATA))
    radds1 = numpy.empty(som1c_1.shape, dtype=numpy.float32)
    radds1[:] = _TARGET_NODATA
    radds1[valid_mask] = (
        rad1p_1_iel[valid_mask] + rad1p_2_iel[valid_mask] *
        ((som1c_1[valid_mask] / som1e_1_iel[valid_mask]) -
            pcemic1_2_iel[valid_mask]))
    rceto2_surface = numpy.empty(som1c_1.shape, dtype=numpy.float32)
    rceto2_surface[:] = _TARGET_NODATA
    rceto2_surface[valid_mask] = numpy.maximum(
        (som1c_1[valid_mask] / som1e_1_iel[valid_mask] + radds1[valid_mask]),
        rad1p_3_iel[valid_mask])
    return rceto2_surface


def calc_tcflow_strucc_1(
        aminrl_1, aminrl_2, strucc_1, struce_1_1, struce_1_2, rnewas_1_1,
        rnewas_2_1, strmax_1, defac, dec1_1, pligst_1, strlig_1, pheff_struc):
    """Calculate total flow out of surface structural C.

    The total potential flow of C out of surface structural material is
    calculated according to its lignin content, the decomposition factor, and
    soil pH. The actual flow is limited by the availability of N and P. N and P
    may be supplied by the mineral source, or by the element (N or P) in the
    decomposing stock.

    Parameters:
        aminrl_1 (numpy.ndarray): derived, average surface mineral N
        aminrl_2 (numpy.ndarray): derived, average surface mineral P
        strucc_1 (numpy.ndarray): state variable, surface structural C
        struce_1_1 (numpy.ndarray): state variable, surface structural N
        struce_1_2 (numpy.ndarray): state variable, surface structural P
        rnewas_1_1 (numpy.ndarray): derived, required C/N ratio for
            aboveground material decomposing to SOM1
        rnewas_2_1 (numpy.ndarray): derived, required C/P ratio for
            aboveground material decomposing to SOM1
        strmax_1 (numpy.ndarray): parameter, maximum decomposition amount
        defac (numpy.ndarray): derived, decomposition factor
        dec1_1 (numpy.ndarray): parameter, maximum decomposition rate
        pligst_1 (numpy.ndarray): parameter, effect of lignin content on
            decomposition rate
        strlig_1 (numpy.ndarray): state variable, lignin content of decomposing
            material
        pheff_struc (numpy.ndarray): derived, effect of soil pH on
            decomposition rate

    Returns:
        tcflow_strucc_1, total flow of C out of surface structural
            material
    """
    valid_mask = (
        (~numpy.isclose(aminrl_1, _SV_NODATA)) &
        (~numpy.isclose(aminrl_2, _SV_NODATA)) &
        (~numpy.isclose(strucc_1, _SV_NODATA)) &
        (~numpy.isclose(struce_1_1, _SV_NODATA)) &
        (~numpy.isclose(struce_1_2, _SV_NODATA)) &
        (rnewas_1_1 != _TARGET_NODATA) &
        (rnewas_2_1 != _TARGET_NODATA) &
        (strmax_1 != _IC_NODATA) &
        (defac != _TARGET_NODATA) &
        (dec1_1 != _IC_NODATA) &
        (pligst_1 != _IC_NODATA) &
        (~numpy.isclose(strlig_1, _SV_NODATA)) &
        (pheff_struc != _TARGET_NODATA))

    potential_flow = numpy.zeros(aminrl_1.shape, dtype=numpy.float32)
    potential_flow[valid_mask] = (
        numpy.minimum(strucc_1[valid_mask], strmax_1[valid_mask]) *
        defac[valid_mask] * dec1_1[valid_mask] *
        numpy.exp(-pligst_1[valid_mask] * strlig_1[valid_mask]) * 0.020833 *
        pheff_struc[valid_mask])

    decompose_mask = (
        ((aminrl_1 > 0.0000001) | ((strucc_1 / struce_1_1) <= rnewas_1_1)) &
        ((aminrl_2 > 0.0000001) | ((strucc_1 / struce_1_2) <= rnewas_2_1)) &
        valid_mask)

    tcflow_strucc_1 = numpy.empty(aminrl_1.shape, dtype=numpy.float32)
    tcflow_strucc_1[:] = _IC_NODATA
    tcflow_strucc_1[valid_mask] = 0.
    tcflow_strucc_1[decompose_mask] = potential_flow[decompose_mask]
    return tcflow_strucc_1


def calc_tcflow_strucc_2(
        aminrl_1, aminrl_2, strucc_2, struce_2_1, struce_2_2, rnewbs_1_1,
        rnewbs_2_1, strmax_2, defac, dec1_2, pligst_2, strlig_2, pheff_struc,
        anerb):
    """Calculate total flow out of soil structural C.

    The total potential flow of C out of soil structural material is
    calculated according to its lignin content, the decomposition factor, and
    soil pH. The actual flow is limited by the availability of N and P. N and P
    may be supplied by the mineral source, or by the element (N or P) in the
    decomposing stock.

    Parameters:
        aminrl_1 (numpy.ndarray): derived, average soil mineral N
        aminrl_2 (numpy.ndarray): derived, average soil mineral P
        strucc_2 (numpy.ndarray): state variable, soil structural C
        struce_2_1 (numpy.ndarray): state variable, soil structural N
        struce_2_2 (numpy.ndarray): state variable, soil structural P
        rnewbs_1_1 (numpy.ndarray): derived, required C/N ratio for
            belowground material decomposing to SOM1
        rnewbs_2_1 (numpy.ndarray): derived, required C/P ratio for
            belowground material decomposing to SOM1
        strmax_2 (numpy.ndarray): parameter, maximum decomposition amount
        defac (numpy.ndarray): derived, decomposition factor
        dec1_2 (numpy.ndarray): parameter, maximum decomposition rate
        pligst_2 (numpy.ndarray): parameter, effect of lignin content on
            decomposition rate
        strlig_2 (numpy.ndarray): state variable, lignin content of decomposing
            material
        pheff_struc (numpy.ndarray): derived, effect of soil pH on
            decomposition rate
        anerb (numpy.ndarray): derived, effect of soil anaerobic conditions on
            decomposition rate

    Returns:
        tcflow_strucc_2, total flow of C out of soil structural
            material
    """
    valid_mask = (
        (~numpy.isclose(aminrl_1, _SV_NODATA)) &
        (~numpy.isclose(aminrl_2, _SV_NODATA)) &
        (~numpy.isclose(strucc_2, _SV_NODATA)) &
        (~numpy.isclose(struce_2_1, _SV_NODATA)) &
        (~numpy.isclose(struce_2_2, _SV_NODATA)) &
        (rnewbs_1_1 != _TARGET_NODATA) &
        (rnewbs_2_1 != _TARGET_NODATA) &
        (strmax_2 != _IC_NODATA) &
        (defac != _TARGET_NODATA) &
        (dec1_2 != _IC_NODATA) &
        (pligst_2 != _IC_NODATA) &
        (~numpy.isclose(strlig_2, _SV_NODATA)) &
        (pheff_struc != _TARGET_NODATA) &
        (anerb != _TARGET_NODATA))

    potential_flow = numpy.zeros(aminrl_1.shape, dtype=numpy.float32)
    potential_flow[valid_mask] = (
        numpy.minimum(strucc_2[valid_mask], strmax_2[valid_mask]) *
        defac[valid_mask] * dec1_2[valid_mask] *
        numpy.exp(-pligst_2[valid_mask] * strlig_2[valid_mask]) * 0.020833 *
        pheff_struc[valid_mask] * anerb[valid_mask])

    decompose_mask = (
        ((aminrl_1 > 0.0000001) | ((strucc_2 / struce_2_1) <= rnewbs_1_1)) &
        ((aminrl_2 > 0.0000001) | ((strucc_2 / struce_2_2) <= rnewbs_2_1)) &
        valid_mask)

    tcflow_strucc_2 = numpy.empty(aminrl_1.shape, dtype=numpy.float32)
    tcflow_strucc_2[:] = _IC_NODATA
    tcflow_strucc_2[valid_mask] = 0.
    tcflow_strucc_2[decompose_mask] = potential_flow[decompose_mask]
    return tcflow_strucc_2


def calc_tcflow_surface(
        aminrl_1, aminrl_2, cstatv, estatv_1, estatv_2, rcetob_1, rcetob_2,
        defac, dec_param, pheff):
    """Calculate total flow of C out of a surface pool.

    The total potential flow of C out of a surface pool is calculated according
    to the decomposition factor and soil pH. The actual flow is limited by the
    availability of N and P. N and P may be supplied by the mineral source, or
    by the element (N or P) in the decomposing stock.

    Parameters:
        aminrl_1 (numpy.ndarray): derived, average surface mineral N
        aminrl_2 (numpy.ndarray): derived, average surface mineral P
        cstatv (numpy.ndarray): state variable, C in decomposing pool
        estatv_1 (numpy.ndarray): state variable, N in decomposing pool
        estatv_2 (numpy.ndarray): state variable, P in decomposing pool
        rcetob_1 (numpy.ndarray): derived, required C/N ratio for
            material entering the receiving pool
        rcetob_2 (numpy.ndarray): derived, required C/P ratio for
            material entering the receiving pool
        defac (numpy.ndarray): derived, decomposition factor
        dec_param (numpy.ndarray): parameter, maximum decomposition rate
        pheff (numpy.ndarray): derived, effect of soil pH on
            decomposition rate

    Returns:
        tcflow, total flow of C out of the decomposing pool
    """
    valid_mask = (
        (~numpy.isclose(aminrl_1, _SV_NODATA)) &
        (~numpy.isclose(aminrl_2, _SV_NODATA)) &
        (~numpy.isclose(cstatv, _SV_NODATA)) &
        (~numpy.isclose(estatv_1, _SV_NODATA)) &
        (~numpy.isclose(estatv_2, _SV_NODATA)) &
        (rcetob_1 != _TARGET_NODATA) &
        (rcetob_2 != _TARGET_NODATA) &
        (defac != _TARGET_NODATA) &
        (dec_param != _IC_NODATA) &
        (pheff != _TARGET_NODATA))

    potential_flow = numpy.zeros(aminrl_1.shape, dtype=numpy.float32)
    potential_flow[valid_mask] = (
        numpy.minimum(
            cstatv[valid_mask] * defac[valid_mask] * dec_param[valid_mask] *
            0.020833 * pheff[valid_mask], cstatv[valid_mask]))

    decompose_mask = (
        ((aminrl_1 > 0.0000001) | ((cstatv / estatv_1) <= rcetob_1)) &
        ((aminrl_2 > 0.0000001) | ((cstatv / estatv_2) <= rcetob_2)) &
        valid_mask)
    tcflow = numpy.empty(aminrl_1.shape, dtype=numpy.float32)
    tcflow[:] = _IC_NODATA
    tcflow[valid_mask] = 0.
    tcflow[decompose_mask] = potential_flow[decompose_mask]
    return tcflow


def calc_tcflow_soil(
        aminrl_1, aminrl_2, cstatv, estatv_1, estatv_2, rcetob_1,
        rcetob_2, defac, dec_param, pheff, anerb):
    """Calculate total flow out of soil metabolic C.

    The total potential flow of C out of soil metabolic material is
    calculated according to the decomposition factor, soil pH, and soil
    anaerobic conditions. The actual flow is limited by the availability of N
    and P. N and P may be supplied by the mineral source, or by the element
    (N or P) in the decomposing stock.

    Parameters:
        aminrl_1 (numpy.ndarray): derived, average soil mineral N
        aminrl_2 (numpy.ndarray): derived, average soil mineral P
        cstatv (numpy.ndarray): state variable, C in decomposing stock
        estatv_1 (numpy.ndarray): state variable, N in decomposing stock
        estatv_2 (numpy.ndarray): state variable, P in decomposing stock
        rcetob_1 (numpy.ndarray): derived, required C/N ratio for
            material entering receiving stock
        rceto1_2 (numpy.ndarray): derived, required C/P ratio for
            material entering receiving stock
        defac (numpy.ndarray): derived, decomposition factor
        dec_param (numpy.ndarray): parameter, maximum decomposition rate
        pheff (numpy.ndarray): derived, effect of soil pH on
            decomposition rate
        anerb (numpy.ndarray): derived, effect of soil anaerobic
            conditions on decomposition rate

    Returns:
        tcflow_soil, total flow of C out of soil metabolic
            material
    """
    valid_mask = (
        (~numpy.isclose(aminrl_1, _SV_NODATA)) &
        (~numpy.isclose(aminrl_2, _SV_NODATA)) &
        (~numpy.isclose(cstatv, _SV_NODATA)) &
        (~numpy.isclose(estatv_1, _SV_NODATA)) &
        (~numpy.isclose(estatv_2, _SV_NODATA)) &
        (rcetob_1 != _TARGET_NODATA) &
        (rcetob_2 != _TARGET_NODATA) &
        (defac != _TARGET_NODATA) &
        (dec_param != _IC_NODATA) &
        (pheff != _TARGET_NODATA) &
        (anerb != _TARGET_NODATA))

    potential_flow = numpy.zeros(aminrl_1.shape, dtype=numpy.float32)
    potential_flow[valid_mask] = (
        numpy.minimum(
            cstatv[valid_mask] * defac[valid_mask] * dec_param[valid_mask] *
            0.020833 * pheff[valid_mask] * anerb[valid_mask],
            cstatv[valid_mask]))

    decompose_mask = (
        ((aminrl_1 > 0.0000001) | ((cstatv / estatv_1) <= rcetob_1)) &
        ((aminrl_2 > 0.0000001) | ((cstatv / estatv_2) <= rcetob_2)) &
        valid_mask)
    tcflow_soil = numpy.empty(aminrl_1.shape, dtype=numpy.float32)
    tcflow_soil[:] = _IC_NODATA
    tcflow_soil[valid_mask] = 0.
    tcflow_soil[decompose_mask] = potential_flow[decompose_mask]
    return tcflow_soil


def calc_tcflow_som1c_2(
        aminrl_1, aminrl_2, som1c_2, som1e_2_1, som1e_2_2, rceto2_1,
        rceto2_2, defac, dec3_2, eftext, anerb, pheff_metab):
    """Calculate total flow out of soil SOM1.

    The total potential flow of C out of soil SOM1 is calculated
    according to the effect of soil texture, anaerobic conditions,
    and soil pH. The actual flow is limited by the availability of N
    and P. N and P may be supplied by the mineral source, or by the
    element (N or P) in the decomposing stock.

    Parameters:
        aminrl_1 (numpy.ndarray): derived, average surface mineral N
        aminrl_2 (numpy.ndarray): derived, average surface mineral P
        som1c_2 (numpy.ndarray): state variable, C in soil SOM1
        som1e_2_1 (numpy.ndarray): state variable, N in soil SOM1
        som1e_2_2 (numpy.ndarray): state variable, P in soil SOM1
        rceto2_1 (numpy.ndarray): derived, required C/N ratio for
            material decomposing to soil SOM2
        rceto2_2 (numpy.ndarray): derived, required C/P ratio for
            material decomposing to soil SOM2
        defac (numpy.ndarray): derived, decomposition factor
        dec3_2 (numpy.ndarray): parameter, maximum decomposition rate
        eftext (numpy.ndarray): derived, effect of soil texture on
            decomposition rate
        anerb (numpy.ndarray): derived, effect of soil anaerobic conditions
            on decomposition rate
        pheff_metab (numpy.ndarray): derived, effect of soil pH on
            decomposition rate

    Returns:
        tcflow_som1c_2, total flow of C out of soil SOM1
    """
    valid_mask = (
        (~numpy.isclose(aminrl_1, _SV_NODATA)) &
        (~numpy.isclose(aminrl_2, _SV_NODATA)) &
        (~numpy.isclose(som1c_2, _SV_NODATA)) &
        (~numpy.isclose(som1e_2_1, _SV_NODATA)) &
        (~numpy.isclose(som1e_2_2, _SV_NODATA)) &
        (rceto2_1 != _TARGET_NODATA) &
        (rceto2_2 != _TARGET_NODATA) &
        (defac != _TARGET_NODATA) &
        (dec3_2 != _IC_NODATA) &
        (eftext != _TARGET_NODATA) &
        (anerb != _TARGET_NODATA) &
        (pheff_metab != _TARGET_NODATA))

    potential_flow = numpy.zeros(aminrl_1.shape, dtype=numpy.float32)
    potential_flow[valid_mask] = (
        som1c_2[valid_mask] * defac[valid_mask] * dec3_2[valid_mask] *
        eftext[valid_mask] * anerb[valid_mask] * 0.020833 *
        pheff_metab[valid_mask])

    decompose_mask = (
        ((aminrl_1 > 0.0000001) | ((som1c_2 / som1e_2_1) <= rceto2_1)) &
        ((aminrl_2 > 0.0000001) | ((som1c_2 / som1e_2_2) <= rceto2_2)) &
        valid_mask)

    tcflow_som1c_2 = numpy.empty(aminrl_1.shape, dtype=numpy.float32)
    tcflow_som1c_2[:] = _IC_NODATA
    tcflow_som1c_2[valid_mask] = 0.
    tcflow_som1c_2[decompose_mask] = potential_flow[decompose_mask]
    return tcflow_som1c_2


def calc_som3_flow(tcflow, fps, animpt, anerb):
    """Calculate the C that flows from soil SOM1 or SOM2 to SOM3.

    The fraction of total flow leaving SOM1 or SOM2 that goes to SOM3 is
    dependent on soil clay content and soil anaerobic conditions.

    Parameters:
        tcflow (numpy.ndarray): derived, total C leaving soil SOM1 or SOM2
        fps (numpy.ndarray): derived, effect of soil clay content on
            decomposition to SOM3
        animpt (numpy.ndarray): parameter, slope of relationship between
            anaerobic conditions and decomposition flow to SOM3
        anerb (numpy.ndarray): derived, impact of soil anaerobic conditions
            on decomposition

    Returns:
        tosom3, C flowing to SOM3
    """
    valid_mask = (
        (tcflow != _IC_NODATA) &
        (fps != _IC_NODATA) &
        (animpt != _IC_NODATA) &
        (anerb != _TARGET_NODATA))
    tosom3 = numpy.empty(tcflow.shape, dtype=numpy.float32)
    tosom3[:] = _IC_NODATA
    tosom3[valid_mask] = (
        tcflow[valid_mask] * fps[valid_mask] *
        (1. + animpt[valid_mask] * (1. - anerb[valid_mask])))
    return tosom3


def calc_som2_flow(som2c_1, cmix, defac):
    """Calculate the C that flows from surface SOM2 to soil SOM2.

    Some C flows from surface SOM2 to soil SOM2 via mixing. This flow is
    controlled by the parameter cmix.

    Parameters:
        som2c_1 (numpy.ndarray): state variable, C in surface SOM2
        cmix (numpy.ndarray): parameter, amount of C flowing via mixing
        defac (numpy.ndarray): derived, decomposition factor

    Returns:
        tcflow, C flowing to soil SOM2 via mixing
    """
    valid_mask = (
        (~numpy.isclose(som2c_1, _SV_NODATA)) &
        (cmix != _IC_NODATA) &
        (defac != _TARGET_NODATA))
    tcflow = numpy.empty(som2c_1.shape, dtype=numpy.float32)
    tcflow[:] = _IC_NODATA
    tcflow[valid_mask] = (
        som2c_1[valid_mask] * cmix[valid_mask] * defac[valid_mask] *
        0.020833)
    return tcflow


def calc_respiration_mineral_flow(cflow, frac_co2, estatv, cstatv):
    """Calculate mineral flow of one element associated with respiration.

    As material decomposes from one stock to another, some CO2 is lost
    to microbial respiration and some nutrient (N or P) moves to the
    mineral pool. Respir.f

    Parameters:
        cflow (numpy.ndarray): derived, C decomposing from one stock
            to another
        frac_co2 (numpy.ndarray): parameter, fraction of decomposing
            C lost as CO2
        estatv (numpy.ndarray): state variable, iel (N or P) in the
            decomposing stock
        cstatv (numpy.ndarray): state variable, C in the decomposing
            stock

    Returns:
        mineral_flow, flow of iel (N or P) accompanying respiration
    """
    valid_mask = (
        (cflow != _IC_NODATA) &
        (frac_co2 != _IC_NODATA) &
        (~numpy.isclose(estatv, _SV_NODATA)) &
        (~numpy.isclose(cstatv, _SV_NODATA)))

    co2_loss = numpy.zeros(cflow.shape, dtype=numpy.float32)
    co2_loss[valid_mask] = cflow[valid_mask] * frac_co2[valid_mask]

    mineral_flow = numpy.empty(cflow.shape, dtype=numpy.float32)
    mineral_flow[:] = _IC_NODATA
    mineral_flow[valid_mask] = (
        co2_loss[valid_mask] * estatv[valid_mask] / cstatv[valid_mask])
    return mineral_flow


def update_gross_mineralization(gross_mineralization, mineral_flow):
    """Update gross N mineralization with current mineral flow.

    Gross mineralization of N during decomposition is used to calculate
    volatilization loss of N after decomposition. It is updated with N
    mineral flow if mineral flow is positive.

    Parameters:
        gross_mineralization (numpy.ndarray): gross N mineralization during
            decomposition
        mineral_flow (numpy.ndarray): N mineral flow

    Returns:
        gromin_updated, updated gross mineralization
    """
    valid_mask = (
        (gross_mineralization != _TARGET_NODATA) &
        (mineral_flow != _IC_NODATA))

    gromin_updated = numpy.empty(
        gross_mineralization.shape, dtype=numpy.float32)
    gromin_updated[:] = _TARGET_NODATA
    gromin_updated[valid_mask] = gross_mineralization[valid_mask]

    update_mask = ((mineral_flow > 0) & valid_mask)
    gromin_updated[update_mask] = (
        gross_mineralization[update_mask] + mineral_flow[update_mask])
    return gromin_updated


def calc_net_cflow(cflow, frac_co2):
    """Calculate net flow of C after loss to CO2.

    As material decomposes from one stock to another, some C is lost to
    CO2 through microbial respiration.  Calculate the net flow of C after
    subtracting losses to CO2.

    Parameters:
        cflow (numpy.ndarray): derived, C decomposing from one stock
            to another
        frac_co2 (numpy.ndarray): parameter, fraction of decomposing
            C lost as CO2

    Returns:
        net_cflow, amount of decomposing C that flows after accounting
            for CO2 losses
    """
    valid_mask = (
        (cflow != _IC_NODATA) &
        (frac_co2 != _IC_NODATA))

    co2_loss = numpy.zeros(cflow.shape, dtype=numpy.float32)
    co2_loss[valid_mask] = cflow[valid_mask] * frac_co2[valid_mask]

    net_cflow = numpy.empty(cflow.shape, dtype=numpy.float32)
    net_cflow[:] = _IC_NODATA
    net_cflow[valid_mask] = cflow[valid_mask] - co2_loss[valid_mask]
    return net_cflow


def calc_net_cflow_tosom2(tcflow, frac_co2, tosom3, cleach):
    """Calculate net flow of C from soil SOM1 to soil SOM2.

    The C flowing from soil SOM1 to SOM2 is the remainder of total flow
    from SOM1, after accounting for losses to CO2 through respiration,
    decomposition to SOM3, and leaching.

    Parameters:
        tcflow (numpy.ndarray): derived, total C decomposing from soil
            SOM1
        frac_co2 (numpy.ndarray): parameter, fraction of decomposing
            C lost as CO2
        tosom3 (numpy.ndarray): derived, C flowing from SOM1 to SOM3
        cleach (numpy.ndarray): derived, leached organic C

    Returns:
        net_tosom2, amount of C that flows from soil SOM1 to soil SOm2
    """
    valid_mask = (
        (tcflow != _IC_NODATA) &
        (frac_co2 != _IC_NODATA) &
        (tosom3 != _IC_NODATA) &
        (cleach != _TARGET_NODATA))
    net_tosom2 = numpy.empty(tcflow.shape, dtype=numpy.float32)
    net_tosom2[:] = _IC_NODATA
    net_tosom2[valid_mask] = (
        tcflow[valid_mask] - (tcflow[valid_mask] * frac_co2[valid_mask]) -
        tosom3[valid_mask] - cleach[valid_mask])
    return net_tosom2


def calc_net_cflow_tosom1(tcflow, frac_co2, tosom3):
    """Calculate net flow of C from soil SOM2 to soil SOM1.

    The C flowing from soil SOM2 to SOM1 is the remainder of total flow
    from SOM2, after accounting for losses to CO2 through respiration
    and decomposition to SOM3.

    Parameters:
        tcflow (numpy.ndarray): derived, total C decomposing from soil
            SOM1
        frac_co2 (numpy.ndarray): parameter, fraction of decomposing
            C lost as CO2
        tosom3 (numpy.ndarray): derived, C flowing from SOM1 to SOM3

    Returns:
        net_tosom1, amount of C that flows from soil SOM2 to soil SOM1
    """
    valid_mask = (
        (tcflow != _IC_NODATA) &
        (frac_co2 != _IC_NODATA) &
        (tosom3 != _IC_NODATA))
    net_tosom1 = numpy.empty(tcflow.shape, dtype=numpy.float32)
    net_tosom1[:] = _IC_NODATA
    net_tosom1[valid_mask] = (
        tcflow[valid_mask] - (tcflow[valid_mask] * frac_co2[valid_mask]) -
        tosom3[valid_mask])
    return net_tosom1


def respiration(
        tcflow_path, frac_co2_path, cstatv_path, estatv_path,
        delta_estatv_path, delta_minerl_1_iel_path, gromin_1_path=None):
    """Calculate and apply flow of N or P during respiration.

    Microbial respiration accompanies decomposition of most stocks.
    Calculate the flow of one element (N or P) to the mineral pool, which
    accompanies this respiration.

    Parameters:
        tcflow_path (string): path to raster containing flow of C that
            is accompanied by respiration
        frac_co2_path (string): path to raster containing fraction of
            C lost to co2
        cstatv_path (string): path to raster containing C state variable
            of decomposing pool
        estatv_path (string): path to raster containing iel (N or P) state
            variable of decomposing pool
        delta_estatv_path (string): path to raster containing change
            in the iel state variable of decomposing pool
        delta_minerl_1_iel_path (string): path to raster containing
            change in surface mineral iel
        gromin_1_path (string): path to raster containing gross
            mineralization of N

    Modifies:
        the raster indicated by `delta_estatv_path`
        the raster indicated by `delta_minerl_1_iel_path`
        the raster indicated by `gromin_1_path`, if supplied

    Returns:
        None
    """
    with tempfile.NamedTemporaryFile(
            prefix='operand_temp', delete=False,
            dir=PROCESSING_DIR) as operand_temp_file:
        operand_temp_path = operand_temp_file.name
    with tempfile.NamedTemporaryFile(
            prefix='d_statv_temp', delete=False,
            dir=PROCESSING_DIR) as d_statv_temp_file:
        d_statv_temp_path = d_statv_temp_file.name

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            tcflow_path, frac_co2_path, estatv_path,
            cstatv_path]],
        calc_respiration_mineral_flow, operand_temp_path, gdal.GDT_Float32,
        _IC_NODATA)
    # mineral flow is removed from the decomposing iel state variable
    shutil.copyfile(delta_estatv_path, d_statv_temp_path)
    raster_difference(
        d_statv_temp_path, _IC_NODATA, operand_temp_path, _IC_NODATA,
        delta_estatv_path, _IC_NODATA)
    # mineral flow is added to surface mineral iel
    shutil.copyfile(delta_minerl_1_iel_path, d_statv_temp_path)
    raster_sum(
        d_statv_temp_path, _IC_NODATA, operand_temp_path, _IC_NODATA,
        delta_minerl_1_iel_path, _IC_NODATA)
    if gromin_1_path:
        shutil.copyfile(gromin_1_path, d_statv_temp_path)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                d_statv_temp_path,
                operand_temp_path]],
            update_gross_mineralization, gromin_1_path,
            gdal.GDT_Float32, _TARGET_NODATA)

    # clean up
    os.remove(operand_temp_path)
    os.remove(d_statv_temp_path)


def nutrient_flow(
        cflow_path, cstatv_donating_path, estatv_donating_path, rcetob_path,
        minerl_1_path, d_estatv_donating_path, d_estatv_receiving_path,
        d_minerl_path, gromin_path=None):
    """Calculate and apply the flow of one nutrient accompanying C.

    As C decomposes from one compartment to another, nutrients (N and P)
    also flow from the donating compartment to the receiving compartment.
    Some N or P may also flow to or from the mineral pool. Calculate and
    apply the flow of iel (N or P) accompanying the given flow of C.

    Parameters:
        cflow_path (string): path to raster containing the flow of C
            from the donating to the receiving pool
        cstatv_donating_path (string): path to raster containing the C
            state variable in the donating pool
        estatv_donating_path (string): path to raster containing the iel
            (N or P) in the donating pool
        rcetob_path (string): path to raster containing required C/iel
            ratio in the receiving pool
        minerl_1_path (string): path to raster containing surface mineral iel
        d_estatv_donating_path (string): path to raster containing change
            in iel in the donating pool
        d_estatv_receiving_path (string): path to raster containing change
            in iel in the receiving pool
        d_minerl_path (string): path to raster containing change in surface
            mineral iel
        gromin_path (string): path to raster containing gross mineralization
            of N

    Modifies:
        the raster indicated by `d_estatv_donating_path`
        the raster indicated by `d_estatv_receiving_path`
        the raster indicated by `d_minerl_path`
        the raster indicated by `gromin_path`, if supplied

    Returns:
        None
    """
    with tempfile.NamedTemporaryFile(
            prefix='operand_temp', delete=False,
            dir=PROCESSING_DIR) as operand_temp_file:
        operand_temp_path = operand_temp_file.name
    with tempfile.NamedTemporaryFile(
            prefix='d_statv_temp', delete=False,
            dir=PROCESSING_DIR) as d_statv_temp_file:
        d_statv_temp_path = d_statv_temp_file.name

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            cflow_path, cstatv_donating_path, rcetob_path,
            estatv_donating_path, minerl_1_path]],
        esched('material_leaving_a'), operand_temp_path, gdal.GDT_Float32,
        _IC_NODATA)
    shutil.copyfile(d_estatv_donating_path, d_statv_temp_path)
    raster_difference(
        d_statv_temp_path, _IC_NODATA, operand_temp_path, _IC_NODATA,
        d_estatv_donating_path, _IC_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            cflow_path, cstatv_donating_path, rcetob_path,
            estatv_donating_path, minerl_1_path]],
        esched('material_arriving_b'), operand_temp_path, gdal.GDT_Float32,
        _IC_NODATA)
    shutil.copyfile(d_estatv_receiving_path, d_statv_temp_path)
    raster_sum(
        d_statv_temp_path, _IC_NODATA, operand_temp_path, _IC_NODATA,
        d_estatv_receiving_path, _IC_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            cflow_path, cstatv_donating_path, rcetob_path,
            estatv_donating_path, minerl_1_path]],
        esched('mineral_flow'), operand_temp_path, gdal.GDT_Float32,
        _IC_NODATA)
    shutil.copyfile(d_minerl_path, d_statv_temp_path)
    raster_sum(
        d_statv_temp_path, _IC_NODATA, operand_temp_path, _IC_NODATA,
        d_minerl_path, _IC_NODATA)
    if gromin_path:
        shutil.copyfile(gromin_path, d_statv_temp_path)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                d_statv_temp_path, operand_temp_path]],
            update_gross_mineralization, gromin_path,
            gdal.GDT_Float32, _TARGET_NODATA)

    # clean up
    os.remove(operand_temp_path)
    os.remove(d_statv_temp_path)


def calc_c_leach(amov_2, tcflow, omlech_3, orglch):
    """Calculate the amount of C leaching from soil SOM1 to stream flow.

    Some C leaches from soil SOM1 if the water flow out of soil layer 2
    is above a critical level.

    Parameters:
        amov_2 (numpy.ndarray): derived, moisture flowing out of soil layer
            2
        tcflow (numpy.ndarray): derived, total flow of C out of soil SOM1
        omlech_3 (numpy.ndarray): parameter, threshold value for amov_2
        orglch (numpy.ndarray): derived, effect of sand content on leaching
            rate

    Returns:
        cleach, C leaching from soil SOM1 to stream flow
    """
    valid_mask = (
        (amov_2 != _TARGET_NODATA) &
        (tcflow != _IC_NODATA) &
        (omlech_3 != _IC_NODATA) &
        (orglch != _IC_NODATA))
    cleach = numpy.empty(amov_2.shape, dtype=numpy.float32)
    cleach[:] = _TARGET_NODATA
    cleach[valid_mask] = 0

    linten = numpy.zeros(amov_2.shape)
    linten[valid_mask] = numpy.minimum(
        (1. - (omlech_3[valid_mask] - amov_2[valid_mask]) /
            omlech_3[valid_mask]), 1.)

    leach_mask = ((amov_2 > 0) & valid_mask)
    cleach[leach_mask] = (
        tcflow[leach_mask] * orglch[leach_mask] * linten[leach_mask])
    return cleach


def remove_leached_iel(
        som1c_2_path, som1e_2_iel_path, cleach_path, d_som1e_2_iel_path,
        iel):
    """Remove N or P leached from soil SOM1.

    As soil SOM1 decomposes into SOM3, some of N and P is lost from SOM1
    through leaching. The amount lost is calculated from the amount of C
    leaching from the soil and the proportion of iel (N or P) in soil SOM1.

    Parameters:
        som1c_2_path (string): path to raster containing C in soil SOM1
        som1e_2_iel_path (string): path to raster containing iel in soil
            SOM1
        cleach_path (string): path to raster containing C leaching from
            SOM1
        d_som1e_2_iel_path (string): path to raster giving change in
            som1e_2_iel
        iel (int): index indicating N (iel == 1) or P (iel == 2))

    Modifies:
        the raster indicated by `d_som1e_2_iel_path`

    Returns:
        None
    """
    def calc_leached_N(som1c_2, som1e_2_1, cleach):
        """Calculate the N leaching from soil SOM1."""
        valid_mask = (
            (~numpy.isclose(som1c_2, _SV_NODATA)) &
            (~numpy.isclose(som1e_2_1, _SV_NODATA)) &
            (cleach != _TARGET_NODATA))
        rceof1_1 = numpy.zeros(som1c_2.shape)
        rceof1_1[valid_mask] = som1c_2[valid_mask] / som1e_2_1[valid_mask] * 2.
        orgflow = numpy.empty(som1c_2.shape, dtype=numpy.float32)
        orgflow[:] = _IC_NODATA
        orgflow[valid_mask] = cleach[valid_mask] / rceof1_1[valid_mask]
        return orgflow

    def calc_leached_P(som1c_2, som1e_2_2, cleach):
        """Calculate the P leaching from soil SOM1."""
        valid_mask = (
            (~numpy.isclose(som1c_2, _SV_NODATA)) &
            (~numpy.isclose(som1e_2_2, _SV_NODATA)) &
            (cleach != _TARGET_NODATA))
        rceof1_2 = numpy.zeros(som1c_2.shape)
        rceof1_2[valid_mask] = (
            som1c_2[valid_mask] / som1e_2_2[valid_mask] * 35.)
        orgflow = numpy.empty(som1c_2.shape, dtype=numpy.float32)
        orgflow[:] = _IC_NODATA
        orgflow[valid_mask] = cleach[valid_mask] / rceof1_2[valid_mask]
        return orgflow

    with tempfile.NamedTemporaryFile(
            prefix='operand_temp', delete=False,
            dir=PROCESSING_DIR) as operand_temp_file:
        operand_temp_path = operand_temp_file.name
    with tempfile.NamedTemporaryFile(
            prefix='d_statv_temp', delete=False,
            dir=PROCESSING_DIR) as d_statv_temp_file:
        d_statv_temp_path = d_statv_temp_file.name

    if iel == 1:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                som1c_2_path, som1e_2_iel_path, cleach_path]],
            calc_leached_N, operand_temp_path,
            gdal.GDT_Float32, _TARGET_NODATA)
    else:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                som1c_2_path, som1e_2_iel_path, cleach_path]],
            calc_leached_P, operand_temp_path,
            gdal.GDT_Float32, _TARGET_NODATA)

    # remove leached iel from SOM1
    shutil.copyfile(d_som1e_2_iel_path, d_statv_temp_path)
    raster_difference(
        d_statv_temp_path, _IC_NODATA,
        operand_temp_path, _IC_NODATA,
        d_som1e_2_iel_path, _IC_NODATA)

    # clean up
    os.remove(operand_temp_path)
    os.remove(d_statv_temp_path)


def calc_pflow(pstatv, rate_param, defac):
    """Calculate the flow of mineral P flowing from one pool to another.

    Mineral P contains multiple pools, including parent material, labile P,
    sorbed and strongly sorbed P, and occluded P. Calculate the flow from one
    mineral P pool to another.

    Parameters:
        pstatv (numpy.ndarray): state variable, P in donating mineral pool
        rate_param (numpy.ndarray): parameter, base rate of flow
        defac (numpy.ndarray): derived, decomposition rate

    Returns:
        pflow, mineral P flowing from donating to receiving pool
    """
    valid_mask = (
        (~numpy.isclose(pstatv, _SV_NODATA)) &
        (rate_param != _IC_NODATA) &
        (defac != _TARGET_NODATA))
    pflow = numpy.empty(pstatv.shape, dtype=numpy.float32)
    pflow[:] = _IC_NODATA
    pflow[valid_mask] = (
        pstatv[valid_mask] * rate_param[valid_mask] * defac[valid_mask] *
        0.020833)
    return pflow


def calc_pflow_to_secndy(minerl_lyr_2, pmnsec_2, fsol, defac):
    """Calculate the flow of mineral to secondary P in one soil layer.

    P flows from the mineral pool of each soil layer into secondary P (strongly
    sorbed P) according to the amount in the mineral pool and the amount of P
    in solution.

    Parameters:
        minerl_lyr_2 (numpy.ndarray): state variable, mineral P in soil layer
            lyr
        pmnsec_2 (numpy.ndarray): parameter, base flow rate
        fsol (numpy.ndarray): derived, fraction of P in solution
        defac (numpy.ndarray): derived, decomposition factor

    Returns:
        fmnsec, flow of mineral P to secondary in one soil layer
    """
    valid_mask = (
        (~numpy.isclose(minerl_lyr_2, _SV_NODATA)) &
        (pmnsec_2 != _IC_NODATA) &
        (fsol != _TARGET_NODATA) &
        (defac != _TARGET_NODATA))
    fmnsec = numpy.empty(minerl_lyr_2.shape, dtype=numpy.float32)
    fmnsec[:] = _IC_NODATA
    fmnsec[valid_mask] = (
        pmnsec_2[valid_mask] * minerl_lyr_2[valid_mask] *
        (1. - fsol[valid_mask]) * defac[valid_mask] * 0.020833)
    return fmnsec


def update_aminrl(
        aminrl_1_path, minerl_1_1_path, aminrl_2_path, minerl_1_2_path,
        fsol_path):
    """Update aminrl_1 and aminrl_2, average mineral N and P in surface soil.

    Aminrl_1, average mineral N, and aminrl_2, average mineral P, represent
    labile N or P available for decomposition. They are kept as a running
    average of the minerl_1_1 (for N) or minerl_1_2 (for P) state variable
    across decomposition time steps.

    Parameters:
        aminrl_1_path (string): path to raster containing average mineral N
        minerl_1_1_path (string): path to raster giving current mineral N
            in soil layer 1
        aminrl_2_path (string): path to raster containing average mineral P
        minerl_1_2_path (string): path to raster giving current mineral N
            in soil layer 2
        fsol_path (string): path to raster giving fraction of mineral P in
            solution

    Modifies:
        the raster indicated by `aminrl_1_path`
        the raster indicated by `aminrl_2_path
    Returns:
        None
    """
    def update_aminrl_1(aminrl_1_prev, minerl_1_1):
        """Update average mineral N."""
        valid_mask = (
            (~numpy.isclose(aminrl_1_prev, _SV_NODATA)) &
            (~numpy.isclose(minerl_1_1, _SV_NODATA)))
        aminrl_1 = numpy.empty(aminrl_1_prev.shape, dtype=numpy.float32)
        aminrl_1[:] = _IC_NODATA
        aminrl_1[valid_mask] = (
            aminrl_1_prev[valid_mask] + minerl_1_1[valid_mask] / 2.)
        return aminrl_1

    def update_aminrl_2(aminrl_2_prev, minerl_1_2, fsol):
        """Update average mineral P.

        Average mineral P is calculated from the fraction of mineral P in
        soil layer 1 that is in solution.

        Parameters:
            aminrl_2_prev (numpy.ndarray): derived, previous average surface
                mineral P
            minerl_1_2 (numpy.ndarray): state variable, current mineral P in
                soil layer 1
            fsol (numpy.ndarray): derived, fraction of labile P in solution

        Returns:
            aminrl_2, updated average mineral P
        """
        valid_mask = (
            (~numpy.isclose(aminrl_2_prev, _SV_NODATA)) &
            (~numpy.isclose(minerl_1_2, _SV_NODATA)) &
            (fsol != _TARGET_NODATA))
        aminrl_2 = numpy.empty(aminrl_2_prev.shape, dtype=numpy.float32)
        aminrl_2[:] = _IC_NODATA
        aminrl_2[valid_mask] = (
            aminrl_2_prev[valid_mask] +
            (minerl_1_2[valid_mask] * fsol[valid_mask]) / 2.)
        return aminrl_2

    with tempfile.NamedTemporaryFile(
            prefix='aminrl_prev', delete=False,
            dir=PROCESSING_DIR) as aminrl_prev_file:
        aminrl_prev_path = aminrl_prev_file.name

    shutil.copyfile(aminrl_1_path, aminrl_prev_path)
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [aminrl_prev_path, minerl_1_1_path]],
        update_aminrl_1, aminrl_1_path, gdal.GDT_Float32, _TARGET_NODATA)

    shutil.copyfile(aminrl_2_path, aminrl_prev_path)
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            aminrl_prev_path, minerl_1_2_path, fsol_path]],
        update_aminrl_2, aminrl_2_path, gdal.GDT_Float32, _TARGET_NODATA)

    # clean up
    os.remove(aminrl_prev_path)


def _decomposition(
        aligned_inputs, current_month, month_index, site_param_table,
        year_reg, month_reg, prev_sv_reg, sv_reg, pp_reg):
    """Update soil C, N and P after decomposition.

    C, N and P move from one surface or soil stock to another depending on the
    availability of N and P in the decomposing stock.  This function covers
    lines 118-323 in Simsom.f, including decomp.f, litdec.f, somdec.f, and
    pschem.f.

    Parameters:
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including precipitation, temperature, and
            site spatial index
        current_month (int): month of the year, such that current_month=1
            indicates January
        month_index (int): month of the simulation, such that month_index=13
            indicates month 13 of the simulation
        site_param_table (dict): map of site spatial indices to dictionaries
            containing site parameters
        year_reg (dict): map of key, path pairs giving paths to annual
            precipitation and annual N deposition rasters
        month_reg (dict): map of key, path pairs giving paths to intermediate
            calculated values that are shared between submodels
        prev_sv_reg (dict): map of key, path pairs giving paths to state
            variables for the previous month
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month
        pp_reg (dict): map of key, path pairs giving persistent parameters
            including required ratios for decomposition, the effect of soil
            texture on decomposition rate, and the effect of soil texture on
            the rate of organic leaching

    Modifies:
        all rasters in sv_reg pertaining to structural, metabolic, som1, som2,
            and som3 C, N, and P; mineral N and P; and parent, secondary, and
            occluded mineral P

    Returns:
        None
    """
    def calc_rprpet(pevap, snowmelt, avh2o_3, precip):
        """Calculate the ratio of precipitation to ref evapotranspiration.

        The ratio of precipitation or snowmelt to reference evapotranspiration
        influences agdefac and bgdefac, the above- and belowground
        decomposition factors.

        Parameters:
            pevap (numpy.ndarray): derived, reference evapotranspiration
            snowmelt (numpy.ndarray): derived, snowmelt occuring this month
            avh2o_3 (numpy.ndarray): derived, moisture in top two soil layers
            precip (numpy.ndarray): input, precipitation for this month

        Returns:
            rprpet, the ratio of precipitation or snowmelt to reference
                evapotranspiration
        """
        valid_mask = (
            (pevap != _TARGET_NODATA) &
            (snowmelt != _TARGET_NODATA) &
            (avh2o_3 != _TARGET_NODATA) &
            (~numpy.isclose(precip, precip_nodata)))

        rprpet = numpy.empty(pevap.shape, dtype=numpy.float32)
        rprpet[:] = _TARGET_NODATA

        snowmelt_mask = (valid_mask & (snowmelt > 0))
        rprpet[snowmelt_mask] = snowmelt[snowmelt_mask] / pevap[snowmelt_mask]

        no_melt_mask = (valid_mask & (snowmelt <= 0))
        rprpet[no_melt_mask] = (
            (avh2o_3[no_melt_mask] + precip[no_melt_mask]) /
            pevap[no_melt_mask])
        return rprpet

    def calc_bgwfunc(rprpet):
        """Calculate the impact of belowground water content on decomposition.

        Bgwfunc reflects the effect of soil moisture on decomposition and is
        also used to calculate shoot senescence due to water stress. It is
        calculated from the ratio of soil water in the top two soil layers to
        reference evapotranspiration.

        Parameters:
            rprpet (numpy.ndarray): derived, ratio of precipitation or snowmelt
                to reference evapotranspiration

        Returns:
            bgwfunc, the effect of soil moisture on decomposition
        """
        valid_mask = (rprpet != _TARGET_NODATA)
        bgwfunc = numpy.empty(rprpet.shape, dtype=numpy.float32)
        bgwfunc[:] = _TARGET_NODATA
        bgwfunc[valid_mask] = (
            1. / (1. + 30 * numpy.exp(-8.5 * rprpet[valid_mask])))
        bgwfunc[(valid_mask & (rprpet > 9))] = 1
        return bgwfunc

    def calc_defac(
            bgwfunc, snow, min_temp, max_temp, teff_1, teff_2, teff_3, teff_4):
        """Calculate decomposition factor.

        The decomposition factor influences the rate of surface and soil
        decomposition and reflects the influence of soil temperature and
        moisture. Lines 151-200, Cycle.f.

        Parameters:
            bgwfunc (numpy.ndarray): derived, effect of soil moisture on
                decomposition
            snow (numpy.ndarray): state variable, current snowpack
            min_temp (numpy.ndarray): input, minimum temperature for the month
            max_temp (numpy.ndarray): input, maximum temperature for the month
            teff_1 (numpy.ndarray): parameter, x location of inflection point
                for calculating the effect of soil temperature on decomposition
                factor
            teff_2 (numpy.ndarray): parameter, y location of inflection point
                for calculating the effect of soil temperature on decomposition
                factor
            teff_3 (numpy.ndarray): parameter, step size for calculating the
                effect of soil temperature on decomposition factor
            teff_4 (numpy.ndarray): parameter, slope of the line at the
                inflection point, for calculating the effect of soil
                temperature on decomposition factor

        Returns:
            defac, aboveground and belowground decomposition factor
        """
        valid_mask = (
            (bgwfunc != _TARGET_NODATA) &
            (snow != _TARGET_NODATA) &
            (~numpy.isclose(min_temp, min_temp_nodata)) &
            (~numpy.isclose(max_temp, max_temp_nodata)) &
            (teff_1 != _IC_NODATA) &
            (teff_2 != _IC_NODATA) &
            (teff_3 != _IC_NODATA) &
            (teff_4 != _IC_NODATA))
        stemp = numpy.empty(bgwfunc.shape, dtype=numpy.float32)
        stemp[:] = _TARGET_NODATA
        stemp[valid_mask] = (min_temp[valid_mask] + max_temp[valid_mask]) / 2.
        stemp[(valid_mask & (snow > 0))] = 0.

        tfunc = numpy.empty(bgwfunc.shape, dtype=numpy.float32)
        tfunc[:] = _TARGET_NODATA
        tfunc[valid_mask] = numpy.maximum(
            0.01,
            (teff_2[valid_mask] + (teff_3[valid_mask] / numpy.pi) *
                numpy.arctan(numpy.pi * teff_4[valid_mask] *
                (stemp[valid_mask] - teff_1[valid_mask]))) /
            (teff_2[valid_mask] + (teff_3[valid_mask] / numpy.pi) *
                numpy.arctan(numpy.pi * teff_4[valid_mask] *
                (30.0 - teff_1[valid_mask]))))

        defac = numpy.empty(bgwfunc.shape, dtype=numpy.float32)
        defac[:] = _TARGET_NODATA
        defac[valid_mask] = numpy.maximum(
            0., tfunc[valid_mask] * bgwfunc[valid_mask])
        return defac

    def calc_pheff_struc(pH):
        """Calculate the effect of soil pH on decomp of structural material.

        The effect of soil pH on decomposition rate is a multiplier ranging
        from 0 to 1.  The effect on decomposition of structural material
        differs from the effect on decomposition of metabolic material in
        the values of two constants.

        Parameters:
            pH (numpy.ndarray): input, soil pH

        Returns:
            pheff_struc, the effect of soil pH on decomposition rate of
                structural material
        """
        valid_mask = (~numpy.isclose(pH, pH_nodata))
        pheff_struc = numpy.empty(pH.shape, dtype=numpy.float32)
        pheff_struc[valid_mask] = numpy.clip(
            (0.5 + (1.1 / numpy.pi) *
                numpy.arctan(numpy.pi * 0.7 * (pH[valid_mask] - 4.))), 0, 1)
        return pheff_struc

    def calc_pheff_metab(pH):
        """Calculate the effect of soil pH on decomp of metabolic material.

        The effect of soil pH on decomposition rate is a multiplier ranging
        from 0 to 1.  The effect on decomposition of structural material
        differs from the effect on decomposition of metabolic material in
        the values of two constants.

        Parameters:
            pH (numpy.ndarray): input, soil pH

        Returns:
            pheff_metab, the effect of soil pH on decomposition rate of
                metabolic material
        """
        valid_mask = (~numpy.isclose(pH, pH_nodata))
        pheff_metab = numpy.empty(pH.shape, dtype=numpy.float32)
        pheff_metab[valid_mask] = numpy.clip(
            (0.5 + (1.14 / numpy.pi) *
                numpy.arctan(numpy.pi * 0.7 * (pH[valid_mask] - 4.8))), 0, 1)
        return pheff_metab

    def calc_pheff_som3(pH):
        """Calculate the effect of soil pH on decomposition of SOM3.

        The effect of soil pH on decomposition rate is a multiplier ranging
        from 0 to 1.  The effect on decomposition of SOM3 differs from the
        effect of pH on decomposition of other pools in the value of
        constants.

        Parameters:
            pH (numpy.ndarray): input, soil pH

        Returns:
            pheff_som3, the effect of soil pH on decomposition rate of
                SOM3
        """
        valid_mask = (~numpy.isclose(pH, pH_nodata))
        pheff_metab = numpy.empty(pH.shape, dtype=numpy.float32)
        pheff_metab[valid_mask] = numpy.clip(
            (0.5 + (1.1 / numpy.pi) *
                numpy.arctan(numpy.pi * 0.7 * (pH[valid_mask] - 3.))), 0, 1)
        return pheff_metab

    precip_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['precip_{}'.format(month_index)])['nodata'][0]
    min_temp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['min_temp_{}'.format(current_month)])['nodata'][0]
    max_temp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['max_temp_{}'.format(current_month)])['nodata'][0]
    pH_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['ph_path'])['nodata'][0]

    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for val in [
            'd_statv_temp', 'operand_temp', 'shwave', 'pevap', 'rprpet',
            'defac', 'anerb', 'gromin_1', 'pheff_struc', 'pheff_metab',
            'aminrl_1', 'aminrl_2', 'fsol', 'tcflow', 'tosom2', 'net_tosom2',
            'tosom1', 'net_tosom1', 'tosom3', 'cleach', 'pheff_som3', 'pflow']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))
        for iel in [1, 2]:
            for val in ['rceto1', 'rceto2', 'rceto3']:
                temp_val_dict['{}_{}'.format(val, iel)] = os.path.join(
                    temp_dir, '{}.tif'.format('{}_{}'.format(val, iel)))
    param_val_dict = {}
    for val in [
            'fwloss_4', 'teff_1', 'teff_2', 'teff_3', 'teff_4', 'drain',
            'aneref_1', 'aneref_2', 'aneref_3', 'sorpmx', 'pslsrb', 'strmax_1',
            'dec1_1', 'pligst_1', 'strmax_2', 'dec1_2', 'pligst_2', 'rsplig',
            'ps1co2_1', 'ps1co2_2', 'dec2_1', 'pcemic1_1_1', 'pcemic1_2_1',
            'pcemic1_3_1', 'pcemic1_1_2', 'pcemic1_2_2', 'pcemic1_3_2',
            'varat1_1_1', 'varat1_2_1', 'varat1_3_1', 'varat1_1_2',
            'varat1_2_2', 'varat1_3_2', 'dec2_2', 'pmco2_1', 'pmco2_2',
            'rad1p_1_1', 'rad1p_2_1', 'rad1p_3_1', 'rad1p_1_2', 'rad1p_2_2',
            'rad1p_3_2', 'dec3_1', 'p1co2a_1', 'varat22_1_1', 'varat22_2_1',
            'varat22_3_1', 'varat22_1_2', 'varat22_2_2', 'varat22_3_2',
            'dec3_2', 'animpt', 'varat3_1_1', 'varat3_2_1', 'varat3_3_1',
            'varat3_1_2', 'varat3_2_2', 'varat3_3_2', 'omlech_3',
            'dec5_2', 'p2co2_2', 'dec5_1', 'p2co2_1', 'dec4', 'p3co2', 'cmix',
            'pparmn_2', 'psecmn_2', 'nlayer', 'pmnsec_2', 'psecoc1', 'psecoc2',
            'vlossg']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (aligned_inputs['site_index'], 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)

    # shwave, shortwave radiation outside the atmosphere
    _shortwave_radiation(
        aligned_inputs['site_index'], current_month, temp_val_dict['shwave'])

    # pet, reference evapotranspiration modified by fwloss parameter
    _reference_evapotranspiration(
        aligned_inputs['max_temp_{}'.format(current_month)],
        aligned_inputs['min_temp_{}'.format(current_month)],
        temp_val_dict['shwave'], param_val_dict['fwloss_4'],
        temp_val_dict['pevap'])

    # rprpet, ratio of precipitation to reference evapotranspiration
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['pevap'], month_reg['snowmelt'],
            sv_reg['avh2o_3_path'],
            aligned_inputs['precip_{}'.format(month_index)]]],
        calc_rprpet, temp_val_dict['rprpet'], gdal.GDT_Float32,
        _TARGET_NODATA)

    # bgwfunc, effect of soil moisture on decomposition
    pygeoprocessing.raster_calculator(
        [(temp_val_dict['rprpet'], 1)],
        calc_bgwfunc, month_reg['bgwfunc'], gdal.GDT_Float32,
        _TARGET_NODATA)

    # defac, decomposition factor calculated from soil temp and moisture
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            month_reg['bgwfunc'], sv_reg['snow_path'],
            aligned_inputs['min_temp_{}'.format(current_month)],
            aligned_inputs['max_temp_{}'.format(current_month)],
            param_val_dict['teff_1'], param_val_dict['teff_2'],
            param_val_dict['teff_3'], param_val_dict['teff_4']]],
        calc_defac, temp_val_dict['defac'], gdal.GDT_Float32,
        _TARGET_NODATA)

    # anerb, impact of soil anaerobic conditions on decomposition
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['rprpet'], temp_val_dict['pevap'],
            param_val_dict['drain'], param_val_dict['aneref_1'],
            param_val_dict['aneref_2'], param_val_dict['aneref_3']]],
        calc_anerb, temp_val_dict['anerb'], gdal.GDT_Float32,
        _TARGET_NODATA)

    # initialize gromin_1, gross mineralization of N
    pygeoprocessing.new_raster_from_base(
        aligned_inputs['site_index'], temp_val_dict['gromin_1'],
        gdal.GDT_Float32, [_TARGET_NODATA], fill_value_list=[0])

    # pH effect on decomposition for structural material
    pygeoprocessing.raster_calculator(
        [(aligned_inputs['ph_path'], 1)],
        calc_pheff_struc, temp_val_dict['pheff_struc'], gdal.GDT_Float32,
        _TARGET_NODATA)

    # pH effect on decomposition for metabolic material
    pygeoprocessing.raster_calculator(
        [(aligned_inputs['ph_path'], 1)],
        calc_pheff_metab, temp_val_dict['pheff_metab'], gdal.GDT_Float32,
        _TARGET_NODATA)

    # initialize aminrl_1 and aminrl_2
    shutil.copyfile(prev_sv_reg['minerl_1_1_path'], temp_val_dict['aminrl_1'])
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            prev_sv_reg['minerl_1_2_path'], param_val_dict['sorpmx'],
            param_val_dict['pslsrb']]],
        fsfunc, temp_val_dict['fsol'], gdal.GDT_Float32, _TARGET_NODATA)
    raster_multiplication(
        prev_sv_reg['minerl_1_2_path'], _SV_NODATA,
        temp_val_dict['fsol'], _TARGET_NODATA,
        temp_val_dict['aminrl_2'], _SV_NODATA)

    _monthly_N_fixation(
        aligned_inputs, month_index, site_param_table, year_reg,
        prev_sv_reg, sv_reg)

    # initialize current month state variables and delta state variable dict
    nlayer_max = int(max(
        site_param_table[i]['nlayer'] for i in site_param_table.iterkeys()))
    delta_sv_dict = {
        'minerl_1_1': os.path.join(temp_dir, 'minerl_1_1.tif'),
        'parent_2': os.path.join(temp_dir, 'parent_2.tif'),
        'secndy_2': os.path.join(temp_dir, 'secndy_2.tif'),
        'occlud': os.path.join(temp_dir, 'occlud.tif'),
    }
    for lyr in xrange(1, nlayer_max + 1):
        state_var = 'minerl_{}_2'.format(lyr)
        shutil.copyfile(
            prev_sv_reg['{}_path'.format(state_var)],
            sv_reg['{}_path'.format(state_var)])
        delta_sv_dict[state_var] = os.path.join(
            temp_dir, '{}.tif'.format(state_var))
    # surface mineral N in sv_reg was initialized by _monthly_N_fixation;
    # copy other layers to current sv_reg here
    for lyr in xrange(2, nlayer_max + 1):
        state_var = 'minerl_{}_1'.format(lyr)
        shutil.copyfile(
            prev_sv_reg['{}_path'.format(state_var)],
            sv_reg['{}_path'.format(state_var)])
    for compartment in ['strlig']:
        for lyr in [1, 2]:
            state_var = '{}_{}'.format(compartment, lyr)
            shutil.copyfile(
                prev_sv_reg['{}_path'.format(state_var)],
                sv_reg['{}_path'.format(state_var)])
    for compartment in ['som3']:
        state_var = '{}c'.format(compartment)
        delta_sv_dict[state_var] = os.path.join(
            temp_dir, '{}.tif'.format(state_var))
        shutil.copyfile(
            prev_sv_reg['{}_path'.format(state_var)],
            sv_reg['{}_path'.format(state_var)])
        for iel in [1, 2]:
            state_var = '{}e_{}'.format(compartment, iel)
            delta_sv_dict[state_var] = os.path.join(
                temp_dir, '{}.tif'.format(state_var))
            shutil.copyfile(
                prev_sv_reg['{}_path'.format(state_var)],
                sv_reg['{}_path'.format(state_var)])
    for compartment in ['struc', 'metab', 'som1', 'som2']:
        for lyr in [1, 2]:
            state_var = '{}c_{}'.format(compartment, lyr)
            delta_sv_dict[state_var] = os.path.join(
                temp_dir, '{}.tif'.format(state_var))
            shutil.copyfile(
                prev_sv_reg['{}_path'.format(state_var)],
                sv_reg['{}_path'.format(state_var)])
            for iel in [1, 2]:
                state_var = '{}e_{}_{}'.format(compartment, lyr, iel)
                delta_sv_dict[state_var] = os.path.join(
                    temp_dir, '{}.tif'.format(state_var))
                shutil.copyfile(
                    prev_sv_reg['{}_path'.format(state_var)],
                    sv_reg['{}_path'.format(state_var)])
    for state_var in ['parent_2', 'secndy_2', 'occlud']:
        shutil.copyfile(
            prev_sv_reg['{}_path'.format(state_var)],
            sv_reg['{}_path'.format(state_var)])

    for _ in xrange(4):
        # initialize change (delta, d) in state variables for this decomp step
        for state_var in delta_sv_dict.iterkeys():
            pygeoprocessing.new_raster_from_base(
                aligned_inputs['site_index'], delta_sv_dict[state_var],
                gdal.GDT_Float32, [_IC_NODATA], fill_value_list=[0])

        # decomposition of structural material in surface and soil
        for lyr in [1, 2]:
            if lyr == 1:
                pygeoprocessing.raster_calculator(
                    [(path, 1) for path in [
                        temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                        sv_reg['strucc_1_path'], sv_reg['struce_1_1_path'],
                        sv_reg['struce_1_2_path'], pp_reg['rnewas_1_1_path'],
                        pp_reg['rnewas_2_1_path'], param_val_dict['strmax_1'],
                        temp_val_dict['defac'], param_val_dict['dec1_1'],
                        param_val_dict['pligst_1'], sv_reg['strlig_1_path'],
                        temp_val_dict['pheff_struc']]],
                    calc_tcflow_strucc_1, temp_val_dict['tcflow'],
                    gdal.GDT_Float32, _IC_NODATA)
            else:
                pygeoprocessing.raster_calculator(
                    [(path, 1) for path in [
                        temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                        sv_reg['strucc_2_path'], sv_reg['struce_2_1_path'],
                        sv_reg['struce_2_2_path'], pp_reg['rnewbs_1_1_path'],
                        pp_reg['rnewbs_2_1_path'], param_val_dict['strmax_2'],
                        temp_val_dict['defac'], param_val_dict['dec1_2'],
                        param_val_dict['pligst_2'], sv_reg['strlig_2_path'],
                        temp_val_dict['pheff_struc'], temp_val_dict['anerb']]],
                    calc_tcflow_strucc_2, temp_val_dict['tcflow'],
                    gdal.GDT_Float32, _IC_NODATA)
            shutil.copyfile(
                delta_sv_dict['strucc_{}'.format(lyr)],
                temp_val_dict['d_statv_temp'])

            raster_difference(
                temp_val_dict['d_statv_temp'], _IC_NODATA,
                temp_val_dict['tcflow'], _IC_NODATA,
                delta_sv_dict['strucc_{}'.format(lyr)], _IC_NODATA)

            # structural material decomposes first to SOM2
            raster_multiplication(
                temp_val_dict['tcflow'], _IC_NODATA,
                sv_reg['strlig_{}_path'.format(lyr)], _SV_NODATA,
                temp_val_dict['tosom2'], _IC_NODATA)
            # microbial respiration with decomposition to SOM2
            respiration(
                temp_val_dict['tosom2'], param_val_dict['rsplig'],
                sv_reg['strucc_{}_path'.format(lyr)],
                sv_reg['struce_{}_1_path'.format(lyr)],
                delta_sv_dict['struce_{}_1'.format(lyr)],
                delta_sv_dict['minerl_1_1'],
                gromin_1_path=temp_val_dict['gromin_1'])
            respiration(
                temp_val_dict['tosom2'], param_val_dict['rsplig'],
                sv_reg['strucc_{}_path'.format(lyr)],
                sv_reg['struce_{}_2_path'.format(lyr)],
                delta_sv_dict['struce_{}_2'.format(lyr)],
                delta_sv_dict['minerl_1_2'])

            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['tosom2'], param_val_dict['rsplig']]],
                calc_net_cflow, temp_val_dict['net_tosom2'], gdal.GDT_Float32,
                _IC_NODATA)
            shutil.copyfile(
                delta_sv_dict['som2c_{}'.format(lyr)],
                temp_val_dict['d_statv_temp'])
            raster_sum(
                temp_val_dict['d_statv_temp'], _IC_NODATA,
                temp_val_dict['net_tosom2'], _IC_NODATA,
                delta_sv_dict['som2c_{}'.format(lyr)], _IC_NODATA)

            if lyr == 1:
                rcetob = 'rnewas'
            else:
                rcetob = 'rnewbs'
            # N and P flows from STRUC to SOM2
            nutrient_flow(
                temp_val_dict['net_tosom2'],
                sv_reg['strucc_{}_path'.format(lyr)],
                sv_reg['struce_{}_1_path'.format(lyr)],
                pp_reg['{}_1_2_path'.format(rcetob)],
                sv_reg['minerl_1_1_path'],
                delta_sv_dict['struce_{}_1'.format(lyr)],
                delta_sv_dict['som2e_{}_1'.format(lyr)],
                delta_sv_dict['minerl_1_1'],
                gromin_path=temp_val_dict['gromin_1'])
            nutrient_flow(
                temp_val_dict['net_tosom2'],
                sv_reg['strucc_{}_path'.format(lyr)],
                sv_reg['struce_{}_2_path'.format(lyr)],
                pp_reg['{}_1_2_path'.format(rcetob)],
                sv_reg['minerl_1_2_path'],
                delta_sv_dict['struce_{}_2'.format(lyr)],
                delta_sv_dict['som2e_{}_2'.format(lyr)],
                delta_sv_dict['minerl_1_2'])

            # structural material decomposes next to SOM1
            raster_difference(
                temp_val_dict['tcflow'], _IC_NODATA, temp_val_dict['tosom2'],
                _IC_NODATA, temp_val_dict['tosom1'], _IC_NODATA)
            # microbial respiration with decomposition to SOM1
            respiration(
                temp_val_dict['tosom1'],
                param_val_dict['ps1co2_{}'.format(lyr)],
                sv_reg['strucc_{}_path'.format(lyr)],
                sv_reg['struce_{}_1_path'.format(lyr)],
                delta_sv_dict['struce_{}_1'.format(lyr)],
                delta_sv_dict['minerl_1_1'],
                gromin_1_path=temp_val_dict['gromin_1'])
            respiration(
                temp_val_dict['tosom1'],
                param_val_dict['ps1co2_{}'.format(lyr)],
                sv_reg['strucc_{}_path'.format(lyr)],
                sv_reg['struce_{}_2_path'.format(lyr)],
                delta_sv_dict['struce_{}_2'.format(lyr)],
                delta_sv_dict['minerl_1_2'])

            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['tosom1'],
                    param_val_dict['ps1co2_{}'.format(lyr)]]],
                calc_net_cflow, temp_val_dict['net_tosom1'], gdal.GDT_Float32,
                _IC_NODATA)
            shutil.copyfile(
                delta_sv_dict['som1c_{}'.format(lyr)],
                temp_val_dict['d_statv_temp'])
            raster_sum(
                temp_val_dict['d_statv_temp'], _IC_NODATA,
                temp_val_dict['net_tosom1'], _IC_NODATA,
                delta_sv_dict['som1c_{}'.format(lyr)], _IC_NODATA)

            if lyr == 1:
                rcetob = 'rnewas'
            else:
                rcetob = 'rnewbs'
            # N and P flows from STRUC to SOM1
            nutrient_flow(
                temp_val_dict['net_tosom1'],
                sv_reg['strucc_{}_path'.format(lyr)],
                sv_reg['struce_{}_1_path'.format(lyr)],
                pp_reg['{}_1_1_path'.format(rcetob)],
                sv_reg['minerl_1_1_path'],
                delta_sv_dict['struce_{}_1'.format(lyr)],
                delta_sv_dict['som1e_{}_1'.format(lyr)],
                delta_sv_dict['minerl_1_1'],
                gromin_path=temp_val_dict['gromin_1'])
            nutrient_flow(
                temp_val_dict['net_tosom1'],
                sv_reg['strucc_{}_path'.format(lyr)],
                sv_reg['struce_{}_2_path'.format(lyr)],
                pp_reg['{}_2_1_path'.format(rcetob)],
                sv_reg['minerl_1_2_path'],
                delta_sv_dict['struce_{}_2'.format(lyr)],
                delta_sv_dict['som1e_{}_2'.format(lyr)],
                delta_sv_dict['minerl_1_2'])

        # decomposition of metabolic material in surface and soil to SOM1
        for lyr in [1, 2]:
            if lyr == 1:
                for iel in [1, 2]:
                    # required ratio for surface metabolic decomposing to SOM1
                    pygeoprocessing.raster_calculator(
                        [(path, 1) for path in [
                            sv_reg['metabe_1_{}_path'.format(iel)],
                            sv_reg['metabc_1_path'],
                            param_val_dict['pcemic1_1_{}'.format(iel)],
                            param_val_dict['pcemic1_2_{}'.format(iel)],
                            param_val_dict['pcemic1_3_{}'.format(iel)]]],
                        _aboveground_ratio,
                        temp_val_dict['rceto1_{}'.format(iel)],
                        gdal.GDT_Float32, _TARGET_NODATA)
                pygeoprocessing.raster_calculator(
                    [(path, 1) for path in [
                        temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                        sv_reg['metabc_1_path'], sv_reg['metabe_1_1_path'],
                        sv_reg['metabe_1_2_path'], temp_val_dict['rceto1_1'],
                        temp_val_dict['rceto1_2'], temp_val_dict['defac'],
                        param_val_dict['dec2_1'],
                        temp_val_dict['pheff_metab']]],
                    calc_tcflow_surface, temp_val_dict['tcflow'],
                    gdal.GDT_Float32, _IC_NODATA)
            else:
                for iel in [1, 2]:
                    # required ratio for soil metabolic decomposing to SOM1
                    pygeoprocessing.raster_calculator(
                        [(path, 1) for path in [
                            temp_val_dict['aminrl_{}'.format(iel)],
                            param_val_dict['varat1_1_{}'.format(iel)],
                            param_val_dict['varat1_2_{}'.format(iel)],
                            param_val_dict['varat1_3_{}'.format(iel)]]],
                        _belowground_ratio,
                        temp_val_dict['rceto1_{}'.format(iel)],
                        gdal.GDT_Float32, _TARGET_NODATA)
                pygeoprocessing.raster_calculator(
                    [(path, 1) for path in [
                        temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                        sv_reg['metabc_2_path'], sv_reg['metabe_2_1_path'],
                        sv_reg['metabe_2_2_path'], temp_val_dict['rceto1_1'],
                        temp_val_dict['rceto1_2'], temp_val_dict['defac'],
                        param_val_dict['dec2_2'], temp_val_dict['pheff_metab'],
                        temp_val_dict['anerb']]],
                    calc_tcflow_soil, temp_val_dict['tcflow'],
                    gdal.GDT_Float32, _IC_NODATA)
            shutil.copyfile(
                delta_sv_dict['metabc_{}'.format(lyr)],
                temp_val_dict['d_statv_temp'])
            raster_difference(
                temp_val_dict['d_statv_temp'], _IC_NODATA,
                temp_val_dict['tcflow'], _IC_NODATA,
                delta_sv_dict['metabc_{}'.format(lyr)], _IC_NODATA)
            # microbial respiration with decomposition to SOM1
            respiration(
                temp_val_dict['tcflow'],
                param_val_dict['pmco2_{}'.format(lyr)],
                sv_reg['metabc_{}_path'.format(lyr)],
                sv_reg['metabe_{}_1_path'.format(lyr)],
                delta_sv_dict['metabe_{}_1'.format(lyr)],
                delta_sv_dict['minerl_1_1'],
                gromin_1_path=temp_val_dict['gromin_1'])
            respiration(
                temp_val_dict['tcflow'],
                param_val_dict['pmco2_{}'.format(lyr)],
                sv_reg['metabc_{}_path'.format(lyr)],
                sv_reg['metabe_{}_2_path'.format(lyr)],
                delta_sv_dict['metabe_{}_2'.format(lyr)],
                delta_sv_dict['minerl_1_2'])

            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['tcflow'],
                    param_val_dict['pmco2_{}'.format(lyr)]]],
                calc_net_cflow, temp_val_dict['net_tosom1'], gdal.GDT_Float32,
                _IC_NODATA)
            shutil.copyfile(
                delta_sv_dict['som1c_{}'.format(lyr)],
                temp_val_dict['d_statv_temp'])
            raster_sum(
                temp_val_dict['d_statv_temp'], _IC_NODATA,
                temp_val_dict['net_tosom1'], _IC_NODATA,
                delta_sv_dict['som1c_{}'.format(lyr)], _IC_NODATA)

            nutrient_flow(
                temp_val_dict['net_tosom1'],
                sv_reg['metabc_{}_path'.format(lyr)],
                sv_reg['metabe_{}_1_path'.format(lyr)],
                temp_val_dict['rceto1_1'], sv_reg['minerl_1_1_path'],
                delta_sv_dict['metabe_{}_1'.format(lyr)],
                delta_sv_dict['som1e_{}_1'.format(lyr)],
                delta_sv_dict['minerl_1_1'],
                gromin_path=temp_val_dict['gromin_1'])
            nutrient_flow(
                temp_val_dict['net_tosom1'],
                sv_reg['metabc_{}_path'.format(lyr)],
                sv_reg['metabe_{}_2_path'.format(lyr)],
                temp_val_dict['rceto1_2'], sv_reg['minerl_1_2_path'],
                delta_sv_dict['metabe_{}_2'.format(lyr)],
                delta_sv_dict['som1e_{}_2'.format(lyr)],
                delta_sv_dict['minerl_1_2'])

        # decomposition of surface SOM1 to surface SOM2: line 63 Somdec.f
        for iel in [1, 2]:
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    sv_reg['som1c_1_path'],
                    sv_reg['som1e_1_{}_path'.format(iel)],
                    param_val_dict['rad1p_1_{}'.format(iel)],
                    param_val_dict['rad1p_2_{}'.format(iel)],
                    param_val_dict['rad1p_3_{}'.format(iel)],
                    param_val_dict['pcemic1_2_{}'.format(iel)]]],
                calc_surface_som2_ratio,
                temp_val_dict['rceto2_{}'.format(iel)],
                gdal.GDT_Float32, _TARGET_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                sv_reg['som1c_1_path'], sv_reg['som1e_1_1_path'],
                sv_reg['som1e_1_2_path'], temp_val_dict['rceto2_1'],
                temp_val_dict['rceto2_2'], temp_val_dict['defac'],
                param_val_dict['dec3_1'],
                temp_val_dict['pheff_struc']]],
            calc_tcflow_surface, temp_val_dict['tcflow'],
            gdal.GDT_Float32, _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som1c_1'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tcflow'], _IC_NODATA,
            delta_sv_dict['som1c_1'], _IC_NODATA)
        # microbial respiration with decomposition to SOM2
        respiration(
            temp_val_dict['tcflow'], param_val_dict['p1co2a_1'],
            sv_reg['som1c_1_path'], sv_reg['som1e_1_1_path'],
            delta_sv_dict['som1e_1_1'], delta_sv_dict['minerl_1_1'],
            gromin_1_path=temp_val_dict['gromin_1'])
        respiration(
            temp_val_dict['tcflow'], param_val_dict['p1co2a_1'],
            sv_reg['som1c_1_path'], sv_reg['som1e_1_2_path'],
            delta_sv_dict['som1e_1_2'], delta_sv_dict['minerl_1_2'])

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['tcflow'], param_val_dict['p1co2a_1']]],
            calc_net_cflow, temp_val_dict['net_tosom2'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som2c_1'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['net_tosom2'], _IC_NODATA,
            delta_sv_dict['som2c_1'], _IC_NODATA)

        # N and P flows from som1e_1 to som2e_1, line 123 Somdec.f
        nutrient_flow(
            temp_val_dict['net_tosom2'], sv_reg['som1c_1_path'],
            sv_reg['som1e_1_1_path'], temp_val_dict['rceto2_1'],
            sv_reg['minerl_1_1_path'], delta_sv_dict['som1e_1_1'],
            delta_sv_dict['som2e_1_1'], delta_sv_dict['minerl_1_1'],
            gromin_path=temp_val_dict['gromin_1'])
        nutrient_flow(
            temp_val_dict['net_tosom2'], sv_reg['som1c_1_path'],
            sv_reg['som1e_1_2_path'], temp_val_dict['rceto2_2'],
            sv_reg['minerl_1_2_path'], delta_sv_dict['som1e_1_2'],
            delta_sv_dict['som2e_1_2'], delta_sv_dict['minerl_1_2'])

        # soil SOM1 decomposes to soil SOM3 and SOM2, line 137 Somdec.f
        for iel in [1, 2]:
            # required ratio for soil SOM1 decomposing to SOM2
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['aminrl_{}'.format(iel)],
                    param_val_dict['varat22_1_{}'.format(iel)],
                    param_val_dict['varat22_2_{}'.format(iel)],
                    param_val_dict['varat22_3_{}'.format(iel)]]],
                _belowground_ratio,
                temp_val_dict['rceto2_{}'.format(iel)],
                gdal.GDT_Float32, _TARGET_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                sv_reg['som1c_2_path'], sv_reg['som1e_2_1_path'],
                sv_reg['som1e_2_2_path'], temp_val_dict['rceto2_1'],
                temp_val_dict['rceto2_2'], temp_val_dict['defac'],
                param_val_dict['dec3_2'], pp_reg['eftext_path'],
                temp_val_dict['anerb'], temp_val_dict['pheff_metab']]],
            calc_tcflow_som1c_2, temp_val_dict['tcflow'],
            gdal.GDT_Float32, _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som1c_2'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tcflow'], _IC_NODATA,
            delta_sv_dict['som1c_2'], _IC_NODATA)
        # microbial respiration with decomposition to SOM3, line 179
        respiration(
            temp_val_dict['tcflow'], pp_reg['p1co2_2_path'],
            sv_reg['som1c_2_path'], sv_reg['som1e_2_1_path'],
            delta_sv_dict['som1e_2_1'], delta_sv_dict['minerl_1_1'],
            gromin_1_path=temp_val_dict['gromin_1'])
        respiration(
            temp_val_dict['tcflow'], pp_reg['p1co2_2_path'],
            sv_reg['som1c_2_path'], sv_reg['som1e_2_2_path'],
            delta_sv_dict['som1e_2_2'], delta_sv_dict['minerl_1_2'])

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['tcflow'], pp_reg['fps1s3_path'],
                param_val_dict['animpt'], temp_val_dict['anerb']]],
            calc_som3_flow, temp_val_dict['tosom3'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(delta_sv_dict['som3c'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tosom3'], _IC_NODATA,
            delta_sv_dict['som3c'], _IC_NODATA)
        for iel in [1, 2]:
            # required ratio for soil SOM1 decomposing to SOM3, line 198
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['aminrl_{}'.format(iel)],
                    param_val_dict['varat3_1_{}'.format(iel)],
                    param_val_dict['varat3_2_{}'.format(iel)],
                    param_val_dict['varat3_3_{}'.format(iel)]]],
                _belowground_ratio,
                temp_val_dict['rceto3_{}'.format(iel)],
                gdal.GDT_Float32, _TARGET_NODATA)
        nutrient_flow(
            temp_val_dict['tosom3'], sv_reg['som1c_2_path'],
            sv_reg['som1e_2_1_path'], temp_val_dict['rceto3_1'],
            sv_reg['minerl_1_1_path'], delta_sv_dict['som1e_2_1'],
            delta_sv_dict['som3e_1'], delta_sv_dict['minerl_1_1'],
            gromin_path=temp_val_dict['gromin_1'])
        nutrient_flow(
            temp_val_dict['tosom3'], sv_reg['som1c_2_path'],
            sv_reg['som1e_2_2_path'], temp_val_dict['rceto3_2'],
            sv_reg['minerl_1_2_path'], delta_sv_dict['som1e_2_2'],
            delta_sv_dict['som3e_2'], delta_sv_dict['minerl_1_2'])

        # organic leaching: line 204 Somdec.f
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                month_reg['amov_2'], temp_val_dict['tcflow'],
                param_val_dict['omlech_3'], pp_reg['orglch_path']]],
            calc_c_leach, temp_val_dict['cleach'], gdal.GDT_Float32,
            _TARGET_NODATA)
        for iel in [1, 2]:
            remove_leached_iel(
                sv_reg['som1c_2_path'], sv_reg['som1e_2_{}_path'.format(iel)],
                temp_val_dict['cleach'],
                delta_sv_dict['som1e_2_{}'.format(iel)], iel)

        # rest of flow from soil SOM1 goes to SOM2
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['tcflow'], pp_reg['p1co2_2_path'],
                temp_val_dict['tosom3'], temp_val_dict['cleach']]],
            calc_net_cflow_tosom2, temp_val_dict['net_tosom2'],
            gdal.GDT_Float32, _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som2c_2'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['net_tosom2'], _IC_NODATA,
            delta_sv_dict['som2c_2'], _IC_NODATA)
        # N and P flows from soil SOM1 to soil SOM2, line 257
        nutrient_flow(
            temp_val_dict['net_tosom2'],
            sv_reg['som1c_2_path'], sv_reg['som1e_2_1_path'],
            temp_val_dict['rceto2_1'], sv_reg['minerl_1_1_path'],
            delta_sv_dict['som1e_2_1'], delta_sv_dict['som2e_2_1'],
            delta_sv_dict['minerl_1_1'],
            gromin_path=temp_val_dict['gromin_1'])
        nutrient_flow(
            temp_val_dict['net_tosom2'],
            sv_reg['som1c_2_path'], sv_reg['som1e_2_2_path'],
            temp_val_dict['rceto2_2'], sv_reg['minerl_1_2_path'],
            delta_sv_dict['som1e_2_2'], delta_sv_dict['som2e_2_2'],
            delta_sv_dict['minerl_1_2'])

        # soil SOM2 decomposing to soil SOM1 and SOM3, line 269
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                sv_reg['som2c_2_path'], sv_reg['som2e_2_1_path'],
                sv_reg['som2e_2_2_path'], temp_val_dict['rceto1_1'],
                temp_val_dict['rceto1_2'], temp_val_dict['defac'],
                param_val_dict['dec5_2'], temp_val_dict['pheff_metab'],
                temp_val_dict['anerb']]],
            calc_tcflow_soil, temp_val_dict['tcflow'],
            gdal.GDT_Float32, _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som2c_2'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tcflow'], _IC_NODATA,
            delta_sv_dict['som2c_2'], _IC_NODATA)
        respiration(
            temp_val_dict['tcflow'], param_val_dict['pmco2_2'],
            sv_reg['som2c_2_path'], sv_reg['som2e_2_1_path'],
            delta_sv_dict['som2e_2_1'], delta_sv_dict['minerl_1_1'],
            gromin_1_path=temp_val_dict['gromin_1'])
        respiration(
            temp_val_dict['tcflow'], param_val_dict['pmco2_2'],
            sv_reg['som2c_2_path'], sv_reg['som2e_2_2_path'],
            delta_sv_dict['som2e_2_2'], delta_sv_dict['minerl_1_2'])

        # soil SOM2 flows first to SOM3
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['tcflow'], pp_reg['fps2s3_path'],
                param_val_dict['animpt'], temp_val_dict['anerb']]],
            calc_som3_flow, temp_val_dict['tosom3'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(delta_sv_dict['som3c'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tosom3'], _IC_NODATA,
            delta_sv_dict['som3c'], _IC_NODATA)
        nutrient_flow(
            temp_val_dict['tosom3'], sv_reg['som2c_2_path'],
            sv_reg['som2e_2_1_path'], temp_val_dict['rceto3_1'],
            sv_reg['minerl_1_1_path'], delta_sv_dict['som2e_2_1'],
            delta_sv_dict['som3e_1'], delta_sv_dict['minerl_1_1'],
            gromin_path=temp_val_dict['gromin_1'])
        nutrient_flow(
            temp_val_dict['tosom3'], sv_reg['som2c_2_path'],
            sv_reg['som2e_2_2_path'], temp_val_dict['rceto3_2'],
            sv_reg['minerl_1_2_path'], delta_sv_dict['som2e_2_2'],
            delta_sv_dict['som3e_2'], delta_sv_dict['minerl_1_2'])

        # rest of flow from soil SOM2 goes to soil SOM1
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['tcflow'], param_val_dict['p2co2_2'],
                temp_val_dict['tosom3']]],
            calc_net_cflow_tosom1, temp_val_dict['net_tosom1'],
            gdal.GDT_Float32, _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som1c_2'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['net_tosom1'], _IC_NODATA,
            delta_sv_dict['som1c_2'], _IC_NODATA)
        nutrient_flow(
            temp_val_dict['net_tosom1'], sv_reg['som2c_2_path'],
            sv_reg['som2e_2_1_path'], temp_val_dict['rceto1_1'],
            sv_reg['minerl_1_1_path'], delta_sv_dict['som2e_2_1'],
            delta_sv_dict['som1e_2_1'], delta_sv_dict['minerl_1_1'],
            gromin_path=temp_val_dict['gromin_1'])
        nutrient_flow(
            temp_val_dict['net_tosom1'], sv_reg['som2c_2_path'],
            sv_reg['som2e_2_2_path'], temp_val_dict['rceto1_2'],
            sv_reg['minerl_1_2_path'], delta_sv_dict['som2e_2_2'],
            delta_sv_dict['som1e_2_2'], delta_sv_dict['minerl_1_2'])

        # surface SOM2 decomposes to surface SOM1
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                sv_reg['som2c_1_path'], sv_reg['som2e_1_1_path'],
                sv_reg['som2e_1_2_path'], temp_val_dict['rceto1_1'],
                temp_val_dict['rceto1_2'], temp_val_dict['defac'],
                param_val_dict['dec5_1'], temp_val_dict['pheff_struc']]],
            calc_tcflow_surface, temp_val_dict['tcflow'],
            gdal.GDT_Float32, _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som2c_1'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tcflow'], _IC_NODATA,
            delta_sv_dict['som2c_1'], _IC_NODATA)
        respiration(
            temp_val_dict['tcflow'], param_val_dict['p2co2_1'],
            sv_reg['som2c_1_path'], sv_reg['som2e_1_1_path'],
            delta_sv_dict['som2e_1_1'], delta_sv_dict['minerl_1_1'],
            gromin_1_path=temp_val_dict['gromin_1'])
        respiration(
            temp_val_dict['tcflow'], param_val_dict['p2co2_1'],
            sv_reg['som2c_1_path'], sv_reg['som2e_1_2_path'],
            delta_sv_dict['som2e_1_2'], delta_sv_dict['minerl_1_2'])

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['tcflow'], param_val_dict['p2co2_1']]],
            calc_net_cflow, temp_val_dict['tosom1'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som1c_1'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tosom1'], _IC_NODATA,
            delta_sv_dict['som1c_1'], _IC_NODATA)
        nutrient_flow(
            temp_val_dict['tosom1'], sv_reg['som2c_1_path'],
            sv_reg['som2e_1_1_path'], temp_val_dict['rceto1_1'],
            sv_reg['minerl_1_1_path'], delta_sv_dict['som2e_1_1'],
            delta_sv_dict['som1e_1_1'], delta_sv_dict['minerl_1_1'],
            gromin_path=temp_val_dict['gromin_1'])
        nutrient_flow(
            temp_val_dict['tosom1'], sv_reg['som2c_1_path'],
            sv_reg['som2e_1_2_path'], temp_val_dict['rceto1_2'],
            sv_reg['minerl_1_2_path'], delta_sv_dict['som2e_1_2'],
            delta_sv_dict['som1e_1_2'], delta_sv_dict['minerl_1_2'])

        # SOM3 decomposing to soil SOM1
        # pH effect on decomposition of SOM3
        pygeoprocessing.raster_calculator(
            [(aligned_inputs['ph_path'], 1)],
            calc_pheff_som3, temp_val_dict['pheff_som3'], gdal.GDT_Float32,
            _TARGET_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['aminrl_1'], temp_val_dict['aminrl_2'],
                sv_reg['som3c_path'], sv_reg['som3e_1_path'],
                sv_reg['som3e_2_path'], temp_val_dict['rceto1_1'],
                temp_val_dict['rceto1_2'], temp_val_dict['defac'],
                param_val_dict['dec4'], temp_val_dict['pheff_som3'],
                temp_val_dict['anerb']]],
            calc_tcflow_soil, temp_val_dict['tcflow'],
            gdal.GDT_Float32, _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som3c'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tcflow'], _IC_NODATA,
            delta_sv_dict['som3c'], _IC_NODATA)
        respiration(
            temp_val_dict['tcflow'], param_val_dict['p3co2'],
            sv_reg['som3c_path'], sv_reg['som3e_1_path'],
            delta_sv_dict['som3e_1'], delta_sv_dict['minerl_1_1'],
            gromin_1_path=temp_val_dict['gromin_1'])
        respiration(
            temp_val_dict['tcflow'], param_val_dict['p3co2'],
            sv_reg['som3c_path'], sv_reg['som3e_2_path'],
            delta_sv_dict['som3e_2'], delta_sv_dict['minerl_1_2'])
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['tcflow'], param_val_dict['p3co2']]],
            calc_net_cflow, temp_val_dict['tosom1'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som1c_2'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tosom1'], _IC_NODATA,
            delta_sv_dict['som1c_2'], _IC_NODATA)

        nutrient_flow(
            temp_val_dict['tosom1'], sv_reg['som3c_path'],
            sv_reg['som3e_1_path'], temp_val_dict['rceto1_1'],
            sv_reg['minerl_1_1_path'], delta_sv_dict['som3e_1'],
            delta_sv_dict['som1e_2_1'], delta_sv_dict['minerl_1_1'],
            gromin_path=temp_val_dict['gromin_1'])
        nutrient_flow(
            temp_val_dict['tosom1'], sv_reg['som3c_path'],
            sv_reg['som3e_2_path'], temp_val_dict['rceto1_2'],
            sv_reg['minerl_1_2_path'], delta_sv_dict['som3e_2'],
            delta_sv_dict['som1e_2_2'], delta_sv_dict['minerl_1_2'])

        # Surface SOM2 flows to soil SOM2 via mixing
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sv_reg['som2c_1_path'], param_val_dict['cmix'],
                temp_val_dict['defac']]],
            calc_som2_flow, temp_val_dict['tcflow'],
            gdal.GDT_Float32, _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som2c_1'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tcflow'], _IC_NODATA,
            delta_sv_dict['som2c_1'], _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['som2c_2'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['tcflow'], _IC_NODATA,
            delta_sv_dict['som2c_2'], _IC_NODATA)
        # ratios for N and P entering soil som2 via mixing
        raster_division(
            sv_reg['som2c_1_path'], _SV_NODATA,
            sv_reg['som2e_1_1_path'], _IC_NODATA,
            temp_val_dict['rceto2_1'], _IC_NODATA)
        raster_division(
            sv_reg['som2c_1_path'], _SV_NODATA,
            sv_reg['som2e_1_2_path'], _IC_NODATA,
            temp_val_dict['rceto2_2'], _IC_NODATA)
        nutrient_flow(
            temp_val_dict['tcflow'], sv_reg['som2c_1_path'],
            sv_reg['som2e_1_1_path'], temp_val_dict['rceto2_1'],
            sv_reg['minerl_1_1_path'], delta_sv_dict['som2e_1_1'],
            delta_sv_dict['som2e_2_1'], delta_sv_dict['minerl_1_1'],
            gromin_path=temp_val_dict['gromin_1'])
        nutrient_flow(
            temp_val_dict['tcflow'], sv_reg['som2c_1_path'],
            sv_reg['som2e_1_2_path'], temp_val_dict['rceto2_2'],
            sv_reg['minerl_1_2_path'], delta_sv_dict['som2e_1_2'],
            delta_sv_dict['som2e_2_2'], delta_sv_dict['minerl_1_2'])

        # P flow from parent to mineral: Pschem.f
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sv_reg['parent_2_path'], param_val_dict['pparmn_2'],
                temp_val_dict['defac']]],
            calc_pflow, temp_val_dict['pflow'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['parent_2'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['pflow'], _IC_NODATA,
            delta_sv_dict['parent_2'], _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['minerl_1_2'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['pflow'], _IC_NODATA,
            delta_sv_dict['minerl_1_2'], _IC_NODATA)

        # P flow from secondary to mineral
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sv_reg['secndy_2_path'], param_val_dict['psecmn_2'],
                temp_val_dict['defac']]],
            calc_pflow, temp_val_dict['pflow'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['secndy_2'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['pflow'], _IC_NODATA,
            delta_sv_dict['secndy_2'], _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['minerl_1_2'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['pflow'], _IC_NODATA,
            delta_sv_dict['minerl_1_2'], _IC_NODATA)

        # P flow from mineral to secondary
        for lyr in xrange(1, nlayer_max + 1):
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    sv_reg['minerl_{}_2_path'.format(lyr)],
                    param_val_dict['pmnsec_2'], temp_val_dict['fsol'],
                    temp_val_dict['defac']]],
                calc_pflow_to_secndy, temp_val_dict['pflow'], gdal.GDT_Float32,
                _IC_NODATA)
            shutil.copyfile(
                delta_sv_dict['minerl_{}_2'.format(lyr)],
                temp_val_dict['d_statv_temp'])
            raster_difference(
                temp_val_dict['d_statv_temp'], _IC_NODATA,
                temp_val_dict['pflow'], _IC_NODATA,
                delta_sv_dict['minerl_{}_2'.format(lyr)], _IC_NODATA)
            shutil.copyfile(
                delta_sv_dict['secndy_2'], temp_val_dict['d_statv_temp'])
            raster_sum(
                temp_val_dict['d_statv_temp'], _IC_NODATA,
                temp_val_dict['pflow'], _IC_NODATA,
                delta_sv_dict['secndy_2'], _IC_NODATA)

        # P flow from secondary to occluded
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sv_reg['secndy_2_path'], param_val_dict['psecoc1'],
                temp_val_dict['defac']]],
            calc_pflow, temp_val_dict['pflow'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['secndy_2'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['pflow'], _IC_NODATA,
            delta_sv_dict['secndy_2'], _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['occlud'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['pflow'], _IC_NODATA,
            delta_sv_dict['occlud'], _IC_NODATA)

        # P flow from occluded to secondary
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sv_reg['occlud_path'], param_val_dict['psecoc2'],
                temp_val_dict['defac']]],
            calc_pflow, temp_val_dict['pflow'], gdal.GDT_Float32,
            _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['occlud'], temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['pflow'], _IC_NODATA,
            delta_sv_dict['occlud'], _IC_NODATA)
        shutil.copyfile(
            delta_sv_dict['secndy_2'], temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _IC_NODATA,
            temp_val_dict['pflow'], _IC_NODATA,
            delta_sv_dict['secndy_2'], _IC_NODATA)

        # accumulate flows
        for compartment in ['struc', 'metab', 'som1', 'som2']:
            for lyr in [1, 2]:
                state_var = '{}c_{}'.format(compartment, lyr)
                shutil.copyfile(
                    sv_reg['{}_path'.format(state_var)],
                    temp_val_dict['operand_temp'])
                raster_sum(
                    delta_sv_dict[state_var], _IC_NODATA,
                    temp_val_dict['operand_temp'], _SV_NODATA,
                    sv_reg['{}_path'.format(state_var)], _SV_NODATA)
                for iel in [1, 2]:
                    state_var = '{}e_{}_{}'.format(compartment, lyr, iel)
                    shutil.copyfile(
                        sv_reg['{}_path'.format(state_var)],
                        temp_val_dict['operand_temp'])
                    raster_sum(
                        delta_sv_dict[state_var], _IC_NODATA,
                        temp_val_dict['operand_temp'], _SV_NODATA,
                        sv_reg['{}_path'.format(state_var)], _SV_NODATA)
        for iel in [1, 2]:
            state_var = 'minerl_1_{}'.format(iel)
            shutil.copyfile(
                sv_reg['{}_path'.format(state_var)],
                temp_val_dict['operand_temp'])
            raster_sum(
                delta_sv_dict[state_var], _IC_NODATA,
                temp_val_dict['operand_temp'], _SV_NODATA,
                sv_reg['{}_path'.format(state_var)], _SV_NODATA)
        for state_var in ['parent_2', 'secndy_2', 'occlud']:
            shutil.copyfile(
                sv_reg['{}_path'.format(state_var)],
                temp_val_dict['operand_temp'])
            raster_sum(
                delta_sv_dict[state_var], _IC_NODATA,
                temp_val_dict['operand_temp'], _SV_NODATA,
                sv_reg['{}_path'.format(state_var)], _SV_NODATA)

        # update aminrl: Simsom.f line 301
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                sv_reg['minerl_1_2_path'], param_val_dict['sorpmx'],
                param_val_dict['pslsrb']]],
            fsfunc, temp_val_dict['fsol'], gdal.GDT_Float32, _TARGET_NODATA)
        update_aminrl(
            temp_val_dict['aminrl_1'], sv_reg['minerl_1_1_path'],
            temp_val_dict['aminrl_2'], sv_reg['minerl_1_2_path'],
            temp_val_dict['fsol'])

    # volatilization loss of N: line 323 Simsom.f
    raster_multiplication(
        temp_val_dict['gromin_1'], _TARGET_NODATA,
        param_val_dict['vlossg'], _IC_NODATA,
        temp_val_dict['operand_temp'], _TARGET_NODATA)
    shutil.copyfile(
        sv_reg['minerl_1_1_path'], temp_val_dict['d_statv_temp'])
    raster_difference(
        temp_val_dict['d_statv_temp'], _SV_NODATA,
        temp_val_dict['operand_temp'], _TARGET_NODATA,
        sv_reg['minerl_1_1_path'], _SV_NODATA)


def partit(
        cpart_path, epart_1_path, epart_2_path, frlign_path, sv_reg,
        site_index_path, site_param_table, lyr):
    """Partition incoming material into structural and metabolic pools.

    When organic material is added to the soil, for example as dead
    biomass falls and becomes litter, or when organic material is added
    from animal waste, it is partitioned into structural (STRUCC_lyr) and
    metabolic (METABC_lyr) material according to the ratio of lignin to N in
    the residue. As residue is partitioned, some N and P may be directly
    absorbed from surface mineral N or P into the residue.

    Parameters:
        cpart_path (string): path to raster containing C in incoming material
            that is to be partitioned
        epart_1_path (string): path to raster containing N in incoming
            material
        epart_2_path (string): path to raster containing P in incoming
            material
        frlign_path (string): path to raster containing fraction of incoming
            material that is lignin
        sv_reg (dict): map of key, path pairs giving paths to current state
            variables
        site_index_path (string): path to site spatial index raster
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        lyr (int): layer which is receiving the incoming material (i.e.,
            1=surface layer, 2=soil layer)

    Modifies:
        the rasters indicated by the following paths:
            sv_reg['minerl_1_1_path']
            sv_reg['minerl_1_2_path']
            sv_reg['metabc_{}_path'.format(lyr)]
            sv_reg['strucc_{}_path'.format(lyr)]
            sv_reg['metabe_{}_1_path'.format(lyr)]
            sv_reg['metabe_{}_2_path'.format(lyr)]
            sv_reg['struce_{}_1_path'.format(lyr)]
            sv_reg['struce_{}_2_path'.format(lyr)]
            sv_reg['strlig_{}_path'.format(lyr)]

    Returns:
        None
    """
    def calc_dirabs(
            cpart, epart_iel, minerl_1_iel, damr_lyr_iel, pabres, damrmn_iel):
        """Calculate direct absorption of mineral N or P.

        When organic material is added to the soil, some mineral N or P may
        be directly absorbed from the surface mineral layer into the incoming
        material. the amount transferred depends on the N or P in the incoming
        material and the required C/N or C/P ratio of receiving material.

        Parameters:
            cpart (numpy.ndarray): derived, C in incoming material
            epart_iel (numpy.ndarray): derived, <iel> in incoming material
            minerl_1_iel (numpy.ndarray): state variable, surface mineral <iel>
            damr_lyr_iel (numpy.ndarray): parameter, fraction of iel in lyr
                absorbed by residue
            pabres (numpy.ndarray): parameter, amount of residue which will
                give maximum direct absorption of iel
            damrmn_iel (numpy.ndarray): parameter, minimum C/iel ratio allowed
                in residue after direct absorption

        Returns:
            dirabs_iel, <iel> (N or P) absorbed from the surface mineral pool
        """
        valid_mask = (
            (cpart != _TARGET_NODATA) &
            (epart_iel != _TARGET_NODATA) &
            (~numpy.isclose(minerl_1_iel, _SV_NODATA)) &
            (damr_lyr_iel != _IC_NODATA) &
            (pabres != _IC_NODATA) &
            (damrmn_iel != _IC_NODATA))

        dirabs_iel = numpy.empty(cpart.shape, dtype=numpy.float32)
        dirabs_iel[:] = _TARGET_NODATA
        dirabs_iel[valid_mask] = 0.

        minerl_mask = ((minerl_1_iel >= 0) & valid_mask)
        dirabs_iel[minerl_mask] = (
            damr_lyr_iel[minerl_mask] * minerl_1_iel[minerl_mask] *
            numpy.maximum(cpart[minerl_mask] / pabres[minerl_mask], 1.))

        # rcetot: C/E ratio of incoming material
        rcetot = numpy.empty(cpart.shape, dtype=numpy.float32)
        rcetot[:] = _IC_NODATA
        e_sufficient_mask = (((epart_iel + dirabs_iel) > 0) & valid_mask)
        rcetot[valid_mask] = 0
        rcetot[e_sufficient_mask] = (
            cpart[e_sufficient_mask] / (
                epart_iel[e_sufficient_mask] + dirabs_iel[e_sufficient_mask]))

        dirabs_mod_mask = ((rcetot < damrmn_iel) & valid_mask)
        dirabs_iel[dirabs_mod_mask] = numpy.maximum(
            cpart[dirabs_mod_mask] / damrmn_iel[dirabs_mod_mask] -
            epart_iel[dirabs_mod_mask], 0.)
        return dirabs_iel

    def calc_d_metabc_lyr(cpart, epart_1, dirabs_1, frlign, spl_1, spl_2):
        """Calculate the change in metabolic C after addition of new material.

        Parameters:
            cpart (numpy.ndarray): C in incoming material
            epart_1 (numpy.ndarray): N in incoming material
            dirabs_1 (numpy.ndarray): derived, direct aborption of mineral N
                into incoming material
            frlign (numpy.ndarray): fraction of incoming material which is
                lignin
            spl_1 (numpy.ndarray): parameter, intercept of regression
                predicting fraction of residue going to metabolic
            spl_2 (numpy.ndarray): parameter, slope of regression predicting
                fraction of residue going to metabolic

        Returns:
            d_metabc_lyr, change in metabc_lyr
        """
        valid_mask = (
            (cpart != _TARGET_NODATA) &
            (epart_1 != _TARGET_NODATA) &
            (dirabs_1 != _TARGET_NODATA) &
            (frlign != _TARGET_NODATA) &
            (spl_1 != _IC_NODATA) &
            (spl_2 != _IC_NODATA))
        movt_mask = ((cpart > 0) & valid_mask)

        # rlnres: ratio of lignin to N in the incoming material
        rlnres = numpy.empty(cpart.shape, dtype=numpy.float32)
        rlnres[:] = _TARGET_NODATA
        rlnres[valid_mask] = 0.
        rlnres[movt_mask] = (
            frlign[movt_mask] / (
                (epart_1[movt_mask] + dirabs_1[movt_mask]) /
                (cpart[movt_mask] * 2.5)))

        # frmet: fraction of cpart that goes to metabolic
        frmet = numpy.empty(cpart.shape, dtype=numpy.float32)
        frmet[:] = _TARGET_NODATA
        frmet[valid_mask] = (
            spl_1[valid_mask] - spl_2[valid_mask] * rlnres[valid_mask])

        lign_exceeded_mask = ((frlign > (1. - frmet)) & valid_mask)
        frmet[lign_exceeded_mask] = 1. - frlign[lign_exceeded_mask]

        d_metabc_lyr = numpy.empty(cpart.shape, dtype=numpy.float32)
        d_metabc_lyr[:] = _TARGET_NODATA
        d_metabc_lyr[valid_mask] = cpart[valid_mask] * frmet[valid_mask]
        return d_metabc_lyr

    def calc_d_strucc_lyr(cpart, d_metabc_lyr):
        """Calculate change in structural C after addition of new material.

        Parameters:
            cpart (numpy.ndarray): derived, C in incoming material
            d_metabc_lyr (numpy.ndarray) derived, change in metabc_lyr

        Returns:
            d_strucc_lyr, change in strucc_lyr
        """
        valid_mask = (
            (cpart != _TARGET_NODATA) &
            (d_metabc_lyr != _TARGET_NODATA))

        d_strucc_lyr = numpy.empty(cpart.shape, dtype=numpy.float32)
        d_strucc_lyr[:] = _TARGET_NODATA
        d_strucc_lyr[valid_mask] = cpart[valid_mask] - d_metabc_lyr[valid_mask]
        return d_strucc_lyr

    def calc_d_struce_lyr_iel(d_strucc_lyr, rcestr_iel):
        """Calculate the change in N or P in structural material in layer lyr.

        Parameters:
            d_strucc_lyr (numpy.ndarray): change in strucc_lyr with addition of
                incoming material
            rcestr_iel (numpy.ndarray): parameter, C/<iel> ratio for structural
                material

        Returns:
            d_struce_lyr_iel, change in structural N or P in layer lyr
        """
        valid_mask = (
            (d_strucc_lyr != _TARGET_NODATA) &
            (rcestr_iel != _IC_NODATA))

        d_struce_lyr_iel = numpy.empty(d_strucc_lyr.shape, dtype=numpy.float32)
        d_struce_lyr_iel[valid_mask] = (
            d_strucc_lyr[valid_mask] / rcestr_iel[valid_mask])
        return d_struce_lyr_iel

    def calc_d_metabe_lyr_iel(cpart, epart_iel, dirabs_iel, d_struce_lyr_iel):
        """Calculate the change in N or P in metabolic material in layer lyr.

        Parameters:
            cpart (numpy.ndarray): C in incoming material
            epart_iel (numpy.ndarray): <iel> in incoming material
            dirabs_iel (numpy.ndarray): <iel> absorbed from the surface mineral
                pool
            d_struce_lyr_iel (numpy.ndarray): change in structural N or P in
                layer lyr

        Returns:
            d_metabe_lyr_iel, change in metabolic N or P in layer lyr
        """
        valid_mask = (
            (cpart != _TARGET_NODATA) &
            (epart_iel != _TARGET_NODATA) &
            (dirabs_iel != _TARGET_NODATA) &
            (d_struce_lyr_iel != _TARGET_NODATA))

        d_metabe_lyr_iel = numpy.empty(cpart.shape, dtype=numpy.float32)
        d_metabe_lyr_iel[:] = _TARGET_NODATA
        d_metabe_lyr_iel[valid_mask] = (
            epart_iel[valid_mask] + dirabs_iel[valid_mask] -
            d_struce_lyr_iel[valid_mask])
        return d_metabe_lyr_iel

    def calc_d_strlig_lyr(frlign, d_strucc_lyr, cpart, strlig_lyr, strucc_lyr):
        """Calculate change in fraction of lignin in structural material.

        Parameters:
            frlign (numpy.ndarray): fraction of incoming material which is
                lignin
            d_strucc_lyr (numpy.ndarray): change in strucc_lyr with addition of
                incoming material
            cpart (numpy.ndarray): C in incoming material
            strlig_lyr (numpy.ndarray): state variable, lignin in structural
                material in receiving layer
            strucc_lyr (numpy.ndarray): state variable, C in structural
                material in layer lyr

        Returns:
            d_strlig_lyr, change in fraction of lignin in structural material
                in layer lyr
        """
        valid_mask = (
            (frlign != _TARGET_NODATA) &
            (d_strucc_lyr != _TARGET_NODATA) &
            (cpart != _TARGET_NODATA) &
            (~numpy.isclose(strlig_lyr, _SV_NODATA)) &
            (~numpy.isclose(strucc_lyr, _SV_NODATA)))
        movt_mask = ((cpart > 0) & valid_mask)

        fligst = numpy.empty(frlign.shape, dtype=numpy.float32)
        fligst[:] = _TARGET_NODATA
        fligst[valid_mask] = 1.
        fligst[movt_mask] = numpy.minimum(
            frlign[movt_mask] / (
                d_strucc_lyr[movt_mask] / cpart[movt_mask]), 1.)

        strlig_lyr_mod = numpy.empty(frlign.shape, dtype=numpy.float32)
        strlig_lyr_mod[:] = _TARGET_NODATA
        strlig_lyr_mod[valid_mask] = (
            ((strlig_lyr[valid_mask] * strucc_lyr[valid_mask]) +
                (fligst[valid_mask] * d_strucc_lyr[valid_mask])) /
            (strucc_lyr[valid_mask] + d_strucc_lyr[valid_mask]))

        d_strlig_lyr = numpy.empty(frlign.shape, dtype=numpy.float32)
        d_strlig_lyr[:] = _IC_NODATA
        d_strlig_lyr[valid_mask] = (
            strlig_lyr_mod[valid_mask] - strlig_lyr[valid_mask])
        return d_strlig_lyr

    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for val in [
            'dirabs_1', 'dirabs_2', 'd_metabc_lyr', 'd_strucc_lyr',
            'd_struce_lyr_iel', 'd_statv_temp', 'operand_temp']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))

    param_val_dict = {}
    for val in [
            'damr_{}_1'.format(lyr), 'damr_{}_2'.format(lyr), 'pabres',
            'damrmn_1', 'damrmn_2', 'spl_1', 'spl_2', 'rcestr_1',
            'rcestr_2']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)

    # direct absorption of N and P from surface mineral layer
    for iel in [1, 2]:
        if iel == 1:
            epart_path = epart_1_path
        else:
            epart_path = epart_2_path
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cpart_path, epart_path,
                sv_reg['minerl_1_{}_path'.format(iel)],
                param_val_dict['damr_{}_{}'.format(lyr, iel)],
                param_val_dict['pabres'],
                param_val_dict['damrmn_{}'.format(iel)]]],
            calc_dirabs, temp_val_dict['dirabs_{}'.format(iel)],
            gdal.GDT_Float32, _TARGET_NODATA)

        # remove direct absorption from surface mineral layer
        shutil.copyfile(
            sv_reg['minerl_1_{}_path'.format(iel)],
            temp_val_dict['d_statv_temp'])
        raster_difference(
            temp_val_dict['d_statv_temp'], _SV_NODATA,
            temp_val_dict['dirabs_{}'.format(iel)], _TARGET_NODATA,
            sv_reg['minerl_1_{}_path'.format(iel)], _SV_NODATA)

    # partition C into structural and metabolic
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            cpart_path, epart_1_path, temp_val_dict['dirabs_1'],
            frlign_path, param_val_dict['spl_1'],
            param_val_dict['spl_2']]],
        calc_d_metabc_lyr, temp_val_dict['d_metabc_lyr'], gdal.GDT_Float32,
        _TARGET_NODATA)
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            cpart_path, temp_val_dict['d_metabc_lyr']]],
        calc_d_strucc_lyr, temp_val_dict['d_strucc_lyr'], gdal.GDT_Float32,
        _TARGET_NODATA)

    shutil.copyfile(
        sv_reg['metabc_{}_path'.format(lyr)], temp_val_dict['d_statv_temp'])
    raster_sum(
        temp_val_dict['d_statv_temp'], _SV_NODATA,
        temp_val_dict['d_metabc_lyr'], _TARGET_NODATA,
        sv_reg['metabc_{}_path'.format(lyr)], _SV_NODATA)
    shutil.copyfile(
        sv_reg['strucc_{}_path'.format(lyr)], temp_val_dict['d_statv_temp'])
    raster_sum(
        temp_val_dict['d_statv_temp'], _SV_NODATA,
        temp_val_dict['d_strucc_lyr'], _TARGET_NODATA,
        sv_reg['strucc_{}_path'.format(lyr)], _SV_NODATA)

    # partition N and P into structural and metabolic
    for iel in [1, 2]:
        if iel == 1:
            epart_path = epart_1_path
        else:
            epart_path = epart_2_path
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['d_strucc_lyr'],
                param_val_dict['rcestr_{}'.format(iel)]]],
            calc_d_struce_lyr_iel, temp_val_dict['d_struce_lyr_iel'],
            gdal.GDT_Float32, _TARGET_NODATA)
        shutil.copyfile(
            sv_reg['struce_{}_{}_path'.format(lyr, iel)],
            temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _SV_NODATA,
            temp_val_dict['d_struce_lyr_iel'], _TARGET_NODATA,
            sv_reg['struce_{}_{}_path'.format(lyr, iel)], _SV_NODATA)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cpart_path, epart_path,
                temp_val_dict['dirabs_{}'.format(iel)],
                temp_val_dict['d_struce_lyr_iel']]],
            calc_d_metabe_lyr_iel, temp_val_dict['operand_temp'],
            gdal.GDT_Float32, _TARGET_NODATA)
        shutil.copyfile(
            sv_reg['metabe_{}_{}_path'.format(lyr, iel)],
            temp_val_dict['d_statv_temp'])
        raster_sum(
            temp_val_dict['d_statv_temp'], _SV_NODATA,
            temp_val_dict['operand_temp'], _TARGET_NODATA,
            sv_reg['metabe_{}_{}_path'.format(lyr, iel)], _SV_NODATA)

    # adjust fraction of lignin in receiving structural pool
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            frlign_path, temp_val_dict['d_strucc_lyr'], cpart_path,
            sv_reg['strlig_{}_path'.format(lyr)],
            sv_reg['strucc_{}_path'.format(lyr)]]],
        calc_d_strlig_lyr, temp_val_dict['operand_temp'], gdal.GDT_Float32,
        _IC_NODATA)
    shutil.copyfile(
        sv_reg['strlig_{}_path'.format(lyr)],
        temp_val_dict['d_statv_temp'])
    raster_sum(
        temp_val_dict['d_statv_temp'], _SV_NODATA,
        temp_val_dict['operand_temp'], _TARGET_NODATA,
        sv_reg['strlig_{}_path'.format(lyr)], _SV_NODATA)


def calc_fall_standing_dead(stdedc, fallrt):
    """Calculate delta C with fall of standing dead.

    Material falls from standing dead biomass into surface litter
    according to a constant monthly fall rate.

    Parameters:
        stdedc (numpy.ndarray): state variable, C in standing dead material
        fallrt (numpy.ndarray): parameter, fraction of standing dead material
            that falls each month

    Returns:
        delta_c_standing_dead, change in C in standing dead
    """
    valid_mask = (
        (~numpy.isclose(stdedc, _SV_NODATA)) &
        (fallrt != _IC_NODATA))
    delta_c_standing_dead = numpy.empty(stdedc.shape, dtype=numpy.float32)
    delta_c_standing_dead[:] = _TARGET_NODATA
    delta_c_standing_dead[valid_mask] = stdedc[valid_mask] * fallrt[valid_mask]
    return delta_c_standing_dead


def calc_root_death(
        average_temperature, rtdtmp, rdr, avh2o_1, deck5, bglivc):
    """Calculate delta C with death of roots.

    Material flows from roots into soil organic matter pools due to root death.
    Root death rate is limited by average temperature and influenced by
    available soil moisture. Change in C is calculated by multiplying the root
    death rate by bglivc, C in live roots.

    Parameters:
        average_temperature (numpy.ndarray): derived, average temperature for
            the current month
        rtdtmp (numpy.ndarray): parameter,  temperature below which root death
            does not occur
        rdr (numpy.ndarray): parameter, maximum root death rate at very dry
            soil conditions
        avh2o_1 (numpy.ndarray): state variable, water available to the current
            plant functional type for growth
        deck5 (numpy.ndarray): parameter, level of available soil water at
            which root death rate is half maximum
        bglivc (numpy.ndarray): state variable, C in belowground live roots

    Returns:
        delta_c_root_death, change in C during root death
    """
    valid_mask = (
        (average_temperature != _IC_NODATA) &
        (rdr != _IC_NODATA) &
        (~numpy.isclose(avh2o_1, _SV_NODATA)) &
        (deck5 != _IC_NODATA) &
        (~numpy.isclose(bglivc, _SV_NODATA)))
    root_death_rate = numpy.empty(bglivc.shape, dtype=numpy.float32)
    root_death_rate[:] = _TARGET_NODATA
    root_death_rate[valid_mask] = 0.

    temp_sufficient_mask = ((average_temperature >= rtdtmp) & valid_mask)
    root_death_rate[temp_sufficient_mask] = numpy.minimum(
        rdr[temp_sufficient_mask] *
        (1.0 - avh2o_1[temp_sufficient_mask] / (
            deck5[temp_sufficient_mask] + avh2o_1[temp_sufficient_mask])),
        0.95)

    delta_c_root_death = numpy.empty(bglivc.shape, dtype=numpy.float32)
    delta_c_root_death[:] = _TARGET_NODATA
    delta_c_root_death[valid_mask] = (
        root_death_rate[valid_mask] * bglivc[valid_mask])
    return delta_c_root_death


def calc_delta_iel(c_state_variable, iel_state_variable, delta_c):
    """Calculate the change in N or P accompanying change in C.

    As C flows out of standing dead biomass or roots, the amount of iel
    (N or P) flowing out of the same pool is calculated from the change in C
    according to the ratio of C to iel in the pool.

    Parameters:
        c_state_variable (numpy.ndarray): state variable, C in the pool that is
            losing material
        iel_state_variable (numpy.ndarray): state variable, N or P in the pool
            that is losing material
        delta_c (numpy.ndarray): derived, change in C. Change in N or P is
            proportional to this amount.

    Returns:
        delta_iel, change in N or P accompanying the change in C
    """
    valid_mask = (
        (~numpy.isclose(c_state_variable, _SV_NODATA)) &
        (~numpy.isclose(iel_state_variable, _SV_NODATA)) &
        (delta_c != _TARGET_NODATA))
    delta_iel = numpy.empty(c_state_variable.shape, dtype=numpy.float32)
    delta_iel[:] = _TARGET_NODATA
    delta_iel[valid_mask] = (
        (iel_state_variable[valid_mask] / c_state_variable[valid_mask]) *
        delta_c[valid_mask])
    return delta_iel


def _death_and_partition(
        state_variable, aligned_inputs, site_param_table, current_month,
        prev_sv_reg, sv_reg, year_reg, pft_id_set, veg_trait_table):
    """Track movement of C, N and P from a pft-level state variable into soil.

    Calculate C, N and P leaving the specified state variable and entering
    surface or soil organic matter pools.  Subtract the change in C, N and P
    from the state variable tracked for each pft, sum up the amounts across
    pfts, and distribute the sum of the material flowing from the state
    variable to surface or soil structural and metabolic pools.

    Parameters:
        state_variable (string): string identifying the state variable that is
            flowing into organic matter.  Must be one of "stded" (for fall
            of standing dead) or "bgliv" (for death of roots). If the state
            variable is stded, material flowing from stded enters surface
            structural and metabolic pools.  If the state variable is bgliv,
            material flowing from bgliv enters soil structural and metabolic
            pools.
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including temperature,
            plant functional type composition, and site spatial index
        site_param_table (dict): map of site spatial indices to dictionaries
            containing site parameters
        current_month (int): month of the year, such that current_month=1
            indicates January
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month
        pft_id_set (set): set of integers identifying plant functional types
        veg_trait_table (dict): map of pft id to dictionaries containing
            plant functional type parameters

    Modifies:
        The rasters indicated by
            sv_reg['<state_variable>c_<pft>_path'] for each pft
            sv_reg['<state_variable>e_1_<pft>_path'] for each pft
            sv_reg['<state_variable>e_2_<pft>_path'] for each pft
            sv_reg['minerl_1_1_path']
            sv_reg['minerl_1_2_path']
            sv_reg['metabc_<lyr>_path']
            sv_reg['strucc_<lyr>_path']
            sv_reg['metabe_<lyr>_1_path']
            sv_reg['metabe_<lyr>_2_path']
            sv_reg['struce_<lyr>_1_path']
            sv_reg['struce_<lyr>_2_path']
            sv_reg['strlig_<lyr>_path']
        where lyr=1 if `state_variable` == 'stded'
              lyr=2 if `state_variable` == 'bgliv'

    Returns:
        None
    """
    def calc_avg_temp(max_temp, min_temp):
        """Calculate average temperature from maximum and minimum temp."""
        valid_mask = (
            (~numpy.isclose(max_temp, max_temp_nodata)) &
            (~numpy.isclose(min_temp, min_temp_nodata)))
        tave = numpy.empty(max_temp.shape, dtype=numpy.float32)
        tave[:] = _IC_NODATA
        tave[valid_mask] = (max_temp[valid_mask] + min_temp[valid_mask]) / 2.
        return tave

    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for val in [
            'tave', 'delta_c', 'delta_iel', 'delta_sv_weighted',
            'operand_temp', 'sum_weighted_delta_C', 'sum_weighted_delta_N',
            'sum_weighted_delta_P', 'weighted_lignin', 'sum_lignin',
            'fraction_lignin']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))
    param_val_dict = {}

    # site-level parameters
    val = 'deck5'
    target_path = os.path.join(temp_dir, '{}.tif'.format(val))
    param_val_dict[val] = target_path
    site_to_val = dict(
        [(site_code, float(table[val])) for
            (site_code, table) in site_param_table.iteritems()])
    pygeoprocessing.reclassify_raster(
        (aligned_inputs['site_index'], 1), site_to_val, target_path,
        gdal.GDT_Float32, _IC_NODATA)

    # pft-level parameters
    for val in['fallrt', 'rtdtmp', 'rdr']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path

    # sum of material across pfts to be partitioned to organic matter
    pygeoprocessing.new_raster_from_base(
        aligned_inputs['site_index'], temp_val_dict['sum_weighted_delta_C'],
        gdal.GDT_Float32, [_TARGET_NODATA], fill_value_list=[0])
    pygeoprocessing.new_raster_from_base(
        aligned_inputs['site_index'], temp_val_dict['sum_weighted_delta_N'],
        gdal.GDT_Float32, [_TARGET_NODATA], fill_value_list=[0])
    pygeoprocessing.new_raster_from_base(
        aligned_inputs['site_index'], temp_val_dict['sum_weighted_delta_P'],
        gdal.GDT_Float32, [_TARGET_NODATA], fill_value_list=[0])
    pygeoprocessing.new_raster_from_base(
        aligned_inputs['site_index'], temp_val_dict['sum_lignin'],
        gdal.GDT_Float32, [_TARGET_NODATA], fill_value_list=[0])

    max_temp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['max_temp_{}'.format(current_month)])['nodata'][0]
    min_temp_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['min_temp_{}'.format(current_month)])['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            aligned_inputs['max_temp_{}'.format(current_month)],
            aligned_inputs['min_temp_{}'.format(current_month)]]],
        calc_avg_temp, temp_val_dict['tave'], gdal.GDT_Float32, _IC_NODATA)

    for pft_i in pft_id_set:
        pft_nodata = pygeoprocessing.get_raster_info(
            aligned_inputs['pft_{}'.format(pft_i)])['nodata'][0]

        # calculate change in C leaving the given state variable
        if state_variable == 'stded':
            fill_val = veg_trait_table[pft_i]['fallrt']
            pygeoprocessing.new_raster_from_base(
                aligned_inputs['site_index'], param_val_dict['fallrt'],
                gdal.GDT_Float32, [_IC_NODATA], fill_value_list=[fill_val])
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    prev_sv_reg['stdedc_{}_path'.format(pft_i)],
                    param_val_dict['fallrt']]],
                calc_fall_standing_dead, temp_val_dict['delta_c'],
                gdal.GDT_Float32, _TARGET_NODATA)
        else:
            for val in ['rtdtmp', 'rdr']:
                fill_val = veg_trait_table[pft_i][val]
                pygeoprocessing.new_raster_from_base(
                    aligned_inputs['site_index'], param_val_dict[val],
                    gdal.GDT_Float32, [_IC_NODATA], fill_value_list=[fill_val])
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['tave'],
                    param_val_dict['rtdtmp'],
                    param_val_dict['rdr'],
                    sv_reg['avh2o_1_{}_path'.format(pft_i)],
                    param_val_dict['deck5'],
                    prev_sv_reg['bglivc_{}_path'.format(pft_i)]]],
                calc_root_death, temp_val_dict['delta_c'],
                gdal.GDT_Float32, _TARGET_NODATA)
        # subtract delta_c from the pft-level state variable
        raster_difference(
            prev_sv_reg['{}c_{}_path'.format(state_variable, pft_i)],
            _SV_NODATA, temp_val_dict['delta_c'], _TARGET_NODATA,
            sv_reg['{}c_{}_path'.format(state_variable, pft_i)], _SV_NODATA)
        # calculate delta C weighted by % cover of this pft
        raster_multiplication(
            temp_val_dict['delta_c'], _TARGET_NODATA,
            aligned_inputs['pft_{}'.format(pft_i)], pft_nodata,
            temp_val_dict['delta_sv_weighted'], _TARGET_NODATA)
        shutil.copyfile(
            temp_val_dict['sum_weighted_delta_C'],
            temp_val_dict['operand_temp'])
        raster_sum(
            temp_val_dict['delta_sv_weighted'], _TARGET_NODATA,
            temp_val_dict['operand_temp'], _TARGET_NODATA,
            temp_val_dict['sum_weighted_delta_C'], _TARGET_NODATA)
        # calculate weighted fraction of flowing C which is lignin
        if state_variable == 'stded':
            frlign_path = year_reg['pltlig_above_{}'.format(pft_i)]
        else:
            frlign_path = year_reg['pltlig_below_{}'.format(pft_i)]
        raster_multiplication(
            temp_val_dict['delta_sv_weighted'], _TARGET_NODATA,
            frlign_path, _TARGET_NODATA,
            temp_val_dict['weighted_lignin'], _TARGET_NODATA)
        shutil.copyfile(
            temp_val_dict['sum_lignin'], temp_val_dict['operand_temp'])
        raster_sum(
            temp_val_dict['weighted_lignin'], _TARGET_NODATA,
            temp_val_dict['operand_temp'], _TARGET_NODATA,
            temp_val_dict['sum_lignin'], _TARGET_NODATA)

        for iel in [1, 2]:
            # calculate N or P flowing out of the pft-level state variable
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    prev_sv_reg['{}c_{}_path'.format(state_variable, pft_i)],
                    prev_sv_reg['{}e_{}_{}_path'.format(
                        state_variable, iel, pft_i)],
                    temp_val_dict['delta_c']]],
                calc_delta_iel, temp_val_dict['delta_iel'],
                gdal.GDT_Float32, _TARGET_NODATA)
            # subtract delta_iel from the pft-level state variable
            raster_difference(
                prev_sv_reg['{}e_{}_{}_path'.format(
                    state_variable, iel, pft_i)], _SV_NODATA,
                temp_val_dict['delta_iel'], _TARGET_NODATA,
                sv_reg['{}e_{}_{}_path'.format(state_variable, iel, pft_i)],
                _SV_NODATA)
            # calculate delta iel weighted by % cover of this pft
            raster_multiplication(
                temp_val_dict['delta_iel'], _TARGET_NODATA,
                aligned_inputs['pft_{}'.format(pft_i)], pft_nodata,
                temp_val_dict['delta_sv_weighted'], _TARGET_NODATA)
            if iel == 1:
                shutil.copyfile(
                    temp_val_dict['sum_weighted_delta_N'],
                    temp_val_dict['operand_temp'])
                raster_sum(
                    temp_val_dict['delta_sv_weighted'], _TARGET_NODATA,
                    temp_val_dict['operand_temp'], _TARGET_NODATA,
                    temp_val_dict['sum_weighted_delta_N'], _TARGET_NODATA)
            else:
                shutil.copyfile(
                    temp_val_dict['sum_weighted_delta_P'],
                    temp_val_dict['operand_temp'])
                raster_sum(
                    temp_val_dict['delta_sv_weighted'], _TARGET_NODATA,
                    temp_val_dict['operand_temp'], _TARGET_NODATA,
                    temp_val_dict['sum_weighted_delta_P'], _TARGET_NODATA)

    # partition sum of C, N and P into structural and metabolic pools
    if state_variable == 'stded':
        lyr = 1
    else:
        lyr = 2
    raster_division(
        temp_val_dict['sum_lignin'], _TARGET_NODATA,
        temp_val_dict['sum_weighted_delta_C'], _TARGET_NODATA,
        temp_val_dict['fraction_lignin'], _TARGET_NODATA)
    partit(
        temp_val_dict['sum_weighted_delta_C'],
        temp_val_dict['sum_weighted_delta_N'],
        temp_val_dict['sum_weighted_delta_P'],
        temp_val_dict['fraction_lignin'],
        sv_reg, aligned_inputs['site_index'], site_param_table, lyr)


def calc_senescence_water_shading(
        aglivc, bgwfunc, fsdeth_1, fsdeth_3, fsdeth_4):
    """Calculate shoot death due to water stress and shading.

    In months where senescence is not scheduled to occur, some shoot death
    may still occur due to water stress and shading.

    Parameters:
        aglivc (numpy.ndarray): state variable, carbon in aboveground live
            biomass
        bgwfunc (numpy.ndarray): derived, effect of soil moisture on
            decomposition and shoot senescence
        fsdeth_1 (numpy.ndarray): parameter, maximum shoot death rate at very
            dry soil conditions
        fsdeth_3 (numpy.ndarray): parameter, additional fraction of shoots
            which die when aglivc is greater than fsdeth_4
        fsdeth_4 (numpy.ndarray): parameter, threshold value for aglivc
            above which shading increases senescence

    Returns:
        fdeth, fraction of aboveground live biomass that is converted to
            standing dead
    """
    valid_mask = (
        (~numpy.isclose(aglivc, _SV_NODATA)) &
        (bgwfunc != _TARGET_NODATA) &
        (fsdeth_1 != _IC_NODATA) &
        (fsdeth_3 != _IC_NODATA) &
        (fsdeth_4 != _IC_NODATA))
    fdeth = numpy.empty(aglivc.shape, dtype=numpy.float32)
    fdeth[:] = _TARGET_NODATA
    fdeth[valid_mask] = fsdeth_1[valid_mask] * (1. - bgwfunc[valid_mask])

    shading_mask = ((aglivc > fsdeth_4) & valid_mask)
    fdeth[shading_mask] = fdeth[shading_mask] + fsdeth_3[shading_mask]
    fdeth[valid_mask] = numpy.minimum(fdeth[valid_mask], 1.)
    return fdeth


def _shoot_senescence(
        pft_id_set, veg_trait_table, prev_sv_reg, sv_reg, month_reg,
        current_month):
    """Senescence of live material to standing dead.

    Live aboveground biomass is converted to standing dead according to
    senescence, which is specified for each pft to occur in one or more months
    of the year. In other months, some senescence may occur because of water
    stress or shading.  During senescence, C, N and P move from agliv to stded
    state variables.

    Parameters:
        pft_id_set (set): set of integers identifying plant functional types
        veg_trait_table (dict): map of pft id to dictionaries containing
            plant functional type parameters
        prev_sv_reg (dict): map of key, path pairs giving paths to state
            variables for the previous month
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month
        month_reg (dict): map of key, path pairs giving paths to intermediate
            calculated values that are shared between submodels
        current_month (int): month of the year, such that current_month=1
            indicates January

    Modifies:
        the rasters indicated by
            sv_reg['aglivc_<pft>_path'] for each pft
            sv_reg['stdedc_<pft>_path'] for each pft
            sv_reg['aglive_1_<pft>_path'] for each pft
            sv_reg['aglive_2_<pft>_path'] for each pft
            sv_reg['crpstg_1_<pft>_path'] for each pft
            sv_reg['crpstg_2_<pft>_path'] for each pft
            sv_reg['stdede_1_<pft>_path'] for each pft
            sv_reg['stdede_2_<pft>_path'] for each pft

    Returns:
        None
    """
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for val in [
            'operand_temp', 'fdeth', 'delta_c', 'delta_iel', 'vol_loss',
            'to_storage', 'to_stdede']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))

    param_val_dict = {}
    for val in[
            'fsdeth_1', 'fsdeth_2', 'fsdeth_3', 'fsdeth_4', 'vlossp',
            'crprtf_1', 'crprtf_2']:
        for pft_i in pft_id_set:
            target_path = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
            param_val_dict['{}_{}'.format(val, pft_i)] = target_path
            fill_val = veg_trait_table[pft_i][val]
            pygeoprocessing.new_raster_from_base(
                prev_sv_reg['aglivc_{}_path'.format(pft_i)], target_path,
                gdal.GDT_Float32, [_IC_NODATA], fill_value_list=[fill_val])

    for pft_i in pft_id_set:
        if str(current_month) == veg_trait_table[pft_i]['senescence_month']:
            temp_val_dict['fdeth'] = param_val_dict[
                'fsdeth_2_{}'.format(pft_i)]
        else:
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    prev_sv_reg['aglivc_{}_path'.format(pft_i)],
                    month_reg['bgwfunc'],
                    param_val_dict['fsdeth_1_{}'.format(pft_i)],
                    param_val_dict['fsdeth_3_{}'.format(pft_i)],
                    param_val_dict['fsdeth_4_{}'.format(pft_i)]]],
                calc_senescence_water_shading, temp_val_dict['fdeth'],
                gdal.GDT_Float32, _TARGET_NODATA)
        # change in C flowing from aboveground live biomass to standing dead
        raster_multiplication(
            temp_val_dict['fdeth'], _TARGET_NODATA,
            prev_sv_reg['aglivc_{}_path'.format(pft_i)], _SV_NODATA,
            temp_val_dict['delta_c'], _TARGET_NODATA)
        raster_difference(
            prev_sv_reg['aglivc_{}_path'.format(pft_i)], _SV_NODATA,
            temp_val_dict['delta_c'], _TARGET_NODATA,
            sv_reg['aglivc_{}_path'.format(pft_i)], _SV_NODATA)
        raster_sum(
            prev_sv_reg['stdedc_{}_path'.format(pft_i)], _SV_NODATA,
            temp_val_dict['delta_c'], _TARGET_NODATA,
            sv_reg['stdedc_{}_path'.format(pft_i)], _SV_NODATA)

        for iel in [1, 2]:
            # change in N or P flowing from aboveground live biomass to dead
            raster_multiplication(
                temp_val_dict['fdeth'], _TARGET_NODATA,
                prev_sv_reg['aglive_{}_{}_path'.format(iel, pft_i)],
                _SV_NODATA, temp_val_dict['delta_iel'], _TARGET_NODATA)
            raster_difference(
                prev_sv_reg['aglive_{}_{}_path'.format(iel, pft_i)],
                _SV_NODATA, temp_val_dict['delta_iel'], _TARGET_NODATA,
                sv_reg['aglive_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)
            if iel == 1:
                # volatilization loss of N
                raster_multiplication(
                    temp_val_dict['delta_iel'], _TARGET_NODATA,
                    param_val_dict['vlossp_{}'.format(pft_i)], _IC_NODATA,
                    temp_val_dict['vol_loss'], _TARGET_NODATA)
                shutil.copyfile(
                    temp_val_dict['delta_iel'], temp_val_dict['operand_temp'])
                raster_difference(
                    temp_val_dict['operand_temp'], _TARGET_NODATA,
                    temp_val_dict['vol_loss'], _TARGET_NODATA,
                    temp_val_dict['delta_iel'], _TARGET_NODATA)
            # a fraction of N and P goes to crop storage
            raster_multiplication(
                temp_val_dict['delta_iel'], _TARGET_NODATA,
                param_val_dict['crprtf_{}_{}'.format(iel, pft_i)], _IC_NODATA,
                temp_val_dict['to_storage'], _TARGET_NODATA)
            raster_sum(
                prev_sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)],
                _SV_NODATA, temp_val_dict['to_storage'], _TARGET_NODATA,
                sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)
            # the rest goes to standing dead biomass
            raster_difference(
                temp_val_dict['delta_iel'], _TARGET_NODATA,
                temp_val_dict['to_storage'], _TARGET_NODATA,
                temp_val_dict['to_stdede'], _TARGET_NODATA)
            raster_sum(
                prev_sv_reg['stdede_{}_{}_path'.format(iel, pft_i)],
                _SV_NODATA, temp_val_dict['to_stdede'], _TARGET_NODATA,
                sv_reg['stdede_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)


def convert_biomass_to_C(biomass_path, c_path):
    """Convert from grams of biomass to grams of carbon.

    The root:shoot submodel calculates potential growth in units of grams of
    biomass, but the growth submodel calculates actual growth from that
    potential growth in units of grams of carbon.  Convert biomass to carbon
    using a conversion factor of 2.5.

    Parameters:
        biomass_path (string): path to raster containing grams of
            biomass
        c_path (string): path to raster that should contain the equivalent
            grams of carbon

    Modifies:
        the raster indicated by `c_path`

    Returns:
        None
    """
    def convert_op(biomass):
        """Convert grams of biomass to grams of carbon."""
        valid_mask = (biomass != _TARGET_NODATA)
        carbon = numpy.empty(biomass.shape, dtype=numpy.float32)
        carbon[:] = _TARGET_NODATA
        carbon[valid_mask] = biomass[valid_mask] / 2.5
        return carbon
    pygeoprocessing.raster_calculator(
        [(biomass_path, 1)], convert_op, c_path, gdal.GDT_Float32,
        _TARGET_NODATA)


def restrict_potential_growth(potenc, availm_1, availm_2, snfxmx_1):
    """Restrict potential growth according to mineral nutrients.

    Limit potential growth by the availability of mineral N and P. Growth only
    occurs if there is some availability of both mineral elements. Line 63
    Restrp.f

    Parameters:
        potenc (numpy.ndarray): potential C production (g C)
        availm_1 (numpy.ndarray): derived, total mineral N available to this
            pft
        availm_2 (numpy.ndarray): derived, total mineral P available to this
            pft
        snfxmx_1 (numpy.ndarray): parameter, maximum symbiotic N fixation rate

    Returns:
        potenc_lim_minerl, potential C production limited by availability of
            mineral nutrients
    """
    valid_mask = (
        (potenc != _TARGET_NODATA) &
        (availm_1 != _TARGET_NODATA) &
        (availm_2 != _TARGET_NODATA) &
        (snfxmx_1 != _IC_NODATA))
    potenc_lim_minerl = numpy.empty(potenc.shape, dtype=numpy.float32)
    potenc_lim_minerl[:] = _TARGET_NODATA
    potenc_lim_minerl[valid_mask] = 0

    growth_mask = (
        ((availm_1 > 0) | (snfxmx_1 > 0)) &
        (availm_2 > 0) &
        valid_mask)
    potenc_lim_minerl[growth_mask] = potenc[growth_mask]
    return potenc_lim_minerl


def c_uptake_aboveground(aglivc, cprodl, rtsh):
    """Do uptake of C from atmosphere to aboveground live biomass.

    Given total C predicted to flow into new growth and the root:shoot ratio
    of new growth, perform the flow of C from the atmosphere into aboveground
    live biomass. Lines 137-146 Growth.f

    Parameters:
        aglivc (numpy.ndarray): state variable, existing C in aboveground live
            biomass
        cprodl (numpy.ndarray): derived, c production limited by nutrient
            availability
        rtsh (numpy.ndarray): derived, root/shoot ratio of new production

    Returns:
        modified_aglivc, modified C in aboveground live biomass
    """
    valid_mask = (
        (~numpy.isclose(aglivc, _SV_NODATA)) &
        (cprodl != _TARGET_NODATA) &
        (rtsh != _TARGET_NODATA))

    c_prod_aboveground = numpy.empty(aglivc.shape, dtype=numpy.float32)
    c_prod_aboveground[:] = _TARGET_NODATA
    c_prod_aboveground[valid_mask] = (
        cprodl[valid_mask] * (1. - rtsh[valid_mask] / (rtsh[valid_mask] + 1.)))

    modified_aglivc = numpy.empty(aglivc.shape, dtype=numpy.float32)
    modified_aglivc[:] = _SV_NODATA
    modified_aglivc[valid_mask] = (
        aglivc[valid_mask] + c_prod_aboveground[valid_mask])
    return modified_aglivc


def c_uptake_belowground(bglivc, cprodl, rtsh):
    """Do uptake of C from atmosphere to belowground live biomass.

    Given total C predicted to flow into new growth and the root:shoot ratio
    of new growth, perform the flow of C from the atmosphere into belowground
    live biomass.  Lines 148-156 Growth.f

    Parameters:
        bglivc (numpy.ndarray): state variable, existing C in belowground live
            biomass
        cprodl (numpy.ndarray): derived, c production limited by nutrient
            availability
        rtsh (numpy.ndarray): derived, root/shoot ratio of new production

    Returns:
        modified_bglivc, modified C in belowground live biomass
    """
    valid_mask = (
        (~numpy.isclose(bglivc, _SV_NODATA)) &
        (cprodl != _TARGET_NODATA) &
        (rtsh != _TARGET_NODATA))

    c_prod_belowground = numpy.empty(bglivc.shape, dtype=numpy.float32)
    c_prod_belowground[:] = _TARGET_NODATA
    c_prod_belowground[valid_mask] = (
        cprodl[valid_mask] * (rtsh[valid_mask] / (rtsh[valid_mask] + 1.)))

    modified_bglivc = numpy.empty(bglivc.shape, dtype=numpy.float32)
    modified_bglivc[:] = _SV_NODATA
    modified_bglivc[valid_mask] = (
        bglivc[valid_mask] + c_prod_belowground[valid_mask])
    return modified_bglivc


def calc_uptake_source(return_type):
    """Calculate uptake of nutrient from available sources."""
    def _uptake(
            eavail_iel, eup_above_iel, eup_below_iel, plantNfix, storage_iel,
            iel):
        """Calculate N or P taken up from one source.

        Given the N or P predicted to flow into new above- and belowground
        production, calculate how much of that nutrient will be taken from the
        crop storage pool and how much will be taken from soil.  For N, some of
        the necessary uptake maybe also come from symbiotic N fixation.

        Parameters:
            eavail_iel (numpy.ndarray): derived, total iel available to this
                plant functional type
            eup_above_iel (numpy.ndarray): derived, iel in new aboveground
                production
            eup_below_iel (numpy.ndarray): derived, iel in new belowground
                production
            plantNfix (numpy.ndarray): derived, symbiotic N fixed by this plant
                functional type
            storage_iel (numpy.ndarray): state variable, iel in crop storage
                pool
            iel (integer): index identifying N or P

        Returns:
            uptake_storage, uptake from crop storage pool, if return_type is
                'uptake_storage'
            uptake_soil, uptake from mineral content of soil layers accessible
                by the plant function type, if return_type is 'uptake_soil'
            uptake_Nfix, uptake from symbiotically fixed nitrogen, if
                return_type is 'uptake_Nfix'
        """
        valid_mask = (
            (eup_above_iel != _TARGET_NODATA) &
            (eup_below_iel != _TARGET_NODATA) &
            (plantNfix != _TARGET_NODATA) &
            (~numpy.isclose(storage_iel, _SV_NODATA)))
        eprodl_iel = numpy.empty(eup_above_iel.shape, dtype=numpy.float32)
        eprodl_iel[:] = _TARGET_NODATA
        eprodl_iel[valid_mask] = (
            eup_above_iel[valid_mask] + eup_below_iel[valid_mask])

        uptake_storage = numpy.empty(eup_above_iel.shape, dtype=numpy.float32)
        uptake_storage[:] = _TARGET_NODATA
        uptake_soil = numpy.empty(eup_above_iel.shape, dtype=numpy.float32)
        uptake_soil[:] = _TARGET_NODATA
        uptake_Nfix = numpy.empty(eup_above_iel.shape, dtype=numpy.float32)
        uptake_Nfix[:] = _TARGET_NODATA

        storage_sufficient_mask = ((eprodl_iel <= storage_iel) & valid_mask)
        uptake_storage[valid_mask] = 0.
        uptake_storage[storage_sufficient_mask] = (
            eprodl_iel[storage_sufficient_mask])
        uptake_soil[storage_sufficient_mask] = 0.
        uptake_Nfix[storage_sufficient_mask] = 0.

        insuff_mask = ((eprodl_iel > storage_iel) & valid_mask)
        uptake_storage[insuff_mask] = storage_iel[insuff_mask]
        if iel == 1:
            uptake_soil[insuff_mask] = numpy.minimum(
                (eprodl_iel[insuff_mask] - storage_iel[insuff_mask] -
                    plantNfix[insuff_mask]),
                (eavail_iel[insuff_mask] - storage_iel[insuff_mask] -
                    plantNfix[insuff_mask]))
            uptake_Nfix[insuff_mask] = plantNfix[insuff_mask]
        else:
            uptake_soil[insuff_mask] = (
                eprodl_iel[insuff_mask] - storage_iel[insuff_mask])

        if return_type == 'uptake_storage':
            return uptake_storage
        elif return_type == 'uptake_soil':
            return uptake_soil
        elif return_type == 'uptake_Nfix':
            return uptake_Nfix
    return _uptake


def calc_aboveground_uptake(total_uptake, eup_above_iel, eup_below_iel):
    """Calculate uptake of nutrient apportioned to aboveground biomass.

    Given the total amount of iel (N or P) taken up from one source, the amount
    of uptake that is apportioned to aboveground biomass is calculated from the
    proportion of demand from above- and belowground growth.

    Parameters:
        total_uptake (numpy.ndarray): derived, uptake of iel from one source
        eup_above_iel (numpy.ndarray): derived, iel in new aboveground growth
        eup_below_iel (numpy.ndarray): derived, iel in new belowground growth

    Returns:
        uptake_above, uptake from one source that is apportioned to aboveground
            biomass
    """
    valid_mask = (
        (total_uptake != _TARGET_NODATA) &
        (eup_above_iel != _TARGET_NODATA) &
        (eup_below_iel != _TARGET_NODATA))
    uptake_above = numpy.empty(total_uptake.shape, dtype=numpy.float32)
    uptake_above[valid_mask] = (
        total_uptake[valid_mask] * (
            eup_above_iel[valid_mask] /
            (eup_above_iel[valid_mask] + eup_below_iel[valid_mask])))
    return uptake_above


def calc_belowground_uptake(total_uptake, eup_above_iel, eup_below_iel):
    """Calculate uptake of nutrient apportioned to _belowground biomass.

    Given the total amount of iel (N or P) taken up from one source, the amount
    of uptake that is apportioned to belowground biomass is calculated from
    the proportion of demand from above- and belowground growth.

    Parameters:
        total_uptake (numpy.ndarray): derived, uptake of iel from one source
        eup_above_iel (numpy.ndarray): derived, iel in new aboveground growth
        eup_below_iel (numpy.ndarray): derived, iel in new belowground growth

    Returns:
        uptake_below, uptake from one source that is apportioned to belowground
            biomass
    """
    valid_mask = (
        (total_uptake != _TARGET_NODATA) &
        (eup_above_iel != _TARGET_NODATA) &
        (eup_below_iel != _TARGET_NODATA))
    uptake_below = numpy.empty(total_uptake.shape, dtype=numpy.float32)
    uptake_below[valid_mask] = (
        total_uptake[valid_mask] * (
            eup_below_iel[valid_mask] /
            (eup_above_iel[valid_mask] + eup_below_iel[valid_mask])))
    return uptake_below


def calc_minerl_uptake_lyr(uptake_soil, minerl_lyr_iel, fsol, availm):
    """Calculate uptake of mineral iel from one soil layer.

    Uptake of mineral iel (N or P) from each soil layer into new growth
    is done according to the proportion of total mineral iel contributed
    by that layer.

    Parameters:
        uptake_soil (numpy.ndarray): derived, total uptake of N or P
            from soil
        minerl_lyr_iel (numpy.ndarray): state variable, mineral iel in
            this soil layer
        fsol (numpy.ndarray): derived, fraction of iel in solution
        availm (numpy.ndarray): derived, sum of mineral iel across
            soil layers accessible by this plant functional type

    Returns:
        minerl_uptake_lyr, uptake of iel from this soil layer
    """
    valid_mask = (
        (uptake_soil != _TARGET_NODATA) &
        (~numpy.isclose(minerl_lyr_iel, _SV_NODATA)) &
        (fsol != _TARGET_NODATA) &
        (availm != _TARGET_NODATA))
    minerl_uptake_lyr = numpy.empty(uptake_soil.shape, dtype=numpy.float32)
    minerl_uptake_lyr[valid_mask] = (
        uptake_soil[valid_mask] * minerl_lyr_iel[valid_mask] *
        fsol[valid_mask] / availm[valid_mask])
    return minerl_uptake_lyr


def nutrient_uptake(
        iel, nlay, percent_cover_path, eup_above_iel_path, eup_below_iel_path,
        plantNfix_path, availm_path, eavail_path, sv_reg, pft_i,
        pslsrb_path, sorpmx_path):
    """Do uptake of N or P from soil and crop storage to aglive and bglive.

    Perform the flows of iel from crop storage pool, soil mineral pools, and
    symbiotic N fixation into new above- and belowground biomass. N and P taken
    up from the soil by one plant functional type are weighted by the percent
    cover of that functional type. Lines 124-156, Restrp.f, lines 186-226,
    Growth.f

    Parameters:
        iel (int): index identifying N (iel=1) or P (iel=2)
        nlay (int): number of soil layers accessible by this plant functional
            type
        percent_cover_path (string): path to raster containing percent cover of
            this plant functional type
        eup_above_iel_path (string): path to raster containing iel in new
            aboveground production
        eup_below_iel_path (string): path to raster containing iel in new
            belowground production
        plantNfix_path (string): path to raster giving symbiotic N fixed by
            this plant functional type
        availm_path (string): path to raster giving the sum of mineral iel
            across soil layers accessible by this plant functional type
        eavail_path (string): path to raster giving total iel available to
            this plant functional type
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month
        pft_i (int): index identifying the current pft
        pslsrb_path (string): path to raster giving pslsrb paramter, slope term
            controlling fraction of mineral P that is labile
        sorpmx_path (string): path to raster giving sorpmx paramter, maximum P
            sorption potential

    Modifies:
        The rasters indicated by
            sv_reg['aglive_<iel>_<pft>_path'], iel in aboveground live biomass,
                for the given iel and current pft
            sv_reg['bglive_<iel>_<pft>_path'], iel in belowground live biomass,
                for the given iel and current pft
            sv_reg['crpstg_<iel>_<pft>_path'], iel in crop storage, for the
                given iel and current pft
            sv_reg['minerl_<layer>_<iel>_path'], iel in the soil mineral pool,
                for each soil layer accessible by the current pft ([1:nlay])

    Returns:
        None
    """
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for val in [
            'uptake_storage', 'uptake_soil', 'uptake_Nfix', 'statv_temp',
            'uptake_above', 'uptake_below', 'fsol', 'minerl_uptake_lyr',
            'uptake_weighted']:
        temp_val_dict[val] = os.path.join(temp_dir, '{}.tif'.format(val))

    # calculate uptake from crop storage
    pft_nodata = pygeoprocessing.get_raster_info(
        percent_cover_path)['nodata'][0]
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            eavail_path, eup_above_iel_path, eup_below_iel_path,
            plantNfix_path, sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)]]] +
        [(iel, 'raw')],
        calc_uptake_source('uptake_storage'), temp_val_dict['uptake_storage'],
        gdal.GDT_Float32, _TARGET_NODATA)
    # calculate uptake from soil
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            eavail_path, eup_above_iel_path, eup_below_iel_path,
            plantNfix_path, sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)]]] +
        [(iel, 'raw')],
        calc_uptake_source('uptake_soil'), temp_val_dict['uptake_soil'],
        gdal.GDT_Float32, _TARGET_NODATA)
    if iel == 1:
        # calculate uptake from symbiotically fixed N
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                eavail_path, eup_above_iel_path, eup_below_iel_path,
                plantNfix_path,
                sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)]]] +
            [(iel, 'raw')],
            calc_uptake_source('uptake_Nfix'), temp_val_dict['uptake_Nfix'],
            gdal.GDT_Float32, _TARGET_NODATA)

    # do uptake from crop storage into aboveground and belowground live
    shutil.copyfile(
        sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)],
        temp_val_dict['statv_temp'])
    raster_difference(
        temp_val_dict['statv_temp'], _SV_NODATA,
        temp_val_dict['uptake_storage'], _TARGET_NODATA,
        sv_reg['crpstg_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['uptake_storage'], eup_above_iel_path,
            eup_below_iel_path]],
        calc_aboveground_uptake, temp_val_dict['uptake_above'],
        gdal.GDT_Float32, _TARGET_NODATA)
    shutil.copyfile(
        sv_reg['aglive_{}_{}_path'.format(iel, pft_i)],
        temp_val_dict['statv_temp'])
    raster_sum(
        temp_val_dict['statv_temp'], _SV_NODATA,
        temp_val_dict['uptake_above'], _TARGET_NODATA,
        sv_reg['aglive_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            temp_val_dict['uptake_storage'], eup_above_iel_path,
            eup_below_iel_path]],
        calc_belowground_uptake, temp_val_dict['uptake_below'],
        gdal.GDT_Float32, _TARGET_NODATA)
    shutil.copyfile(
        sv_reg['bglive_{}_{}_path'.format(iel, pft_i)],
        temp_val_dict['statv_temp'])
    raster_sum(
        temp_val_dict['statv_temp'], _SV_NODATA,
        temp_val_dict['uptake_below'], _TARGET_NODATA,
        sv_reg['bglive_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)

    # uptake from each soil layer in proportion to its contribution to availm
    for lyr in xrange(1, nlay + 1):
        if iel == 2:
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    sv_reg['minerl_1_2_path'], sorpmx_path,
                    pslsrb_path]],
                fsfunc, temp_val_dict['fsol'], gdal.GDT_Float32,
                _TARGET_NODATA)
        else:
            pygeoprocessing.new_raster_from_base(
                sv_reg['aglive_{}_{}_path'.format(iel, pft_i)],
                temp_val_dict['fsol'],
                gdal.GDT_Float32, [_IC_NODATA], fill_value_list=[1.])
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['uptake_soil'],
                sv_reg['minerl_{}_{}_path'.format(lyr, iel)],
                temp_val_dict['fsol'], availm_path]],
            calc_minerl_uptake_lyr, temp_val_dict['minerl_uptake_lyr'],
            gdal.GDT_Float32, _TARGET_NODATA)

        # uptake removed from soil is weighted by pft % cover
        raster_multiplication(
            percent_cover_path, pft_nodata,
            temp_val_dict['minerl_uptake_lyr'], _TARGET_NODATA,
            temp_val_dict['uptake_weighted'], _TARGET_NODATA)
        shutil.copyfile(
            sv_reg['minerl_{}_{}_path'.format(lyr, iel)],
            temp_val_dict['statv_temp'])
        raster_difference(
            temp_val_dict['statv_temp'], _SV_NODATA,
            temp_val_dict['uptake_weighted'], _TARGET_NODATA,
            sv_reg['minerl_{}_{}_path'.format(lyr, iel)], _SV_NODATA)

        # uptake from minerl iel in lyr into above and belowground live
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['minerl_uptake_lyr'], eup_above_iel_path,
                eup_below_iel_path]],
            calc_aboveground_uptake, temp_val_dict['uptake_above'],
            gdal.GDT_Float32, _TARGET_NODATA)
        shutil.copyfile(
            sv_reg['aglive_{}_{}_path'.format(iel, pft_i)],
            temp_val_dict['statv_temp'])
        raster_sum(
            temp_val_dict['statv_temp'], _SV_NODATA,
            temp_val_dict['uptake_above'], _TARGET_NODATA,
            sv_reg['aglive_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['minerl_uptake_lyr'], eup_above_iel_path,
                eup_below_iel_path]],
            calc_belowground_uptake, temp_val_dict['uptake_below'],
            gdal.GDT_Float32, _TARGET_NODATA)
        shutil.copyfile(
            sv_reg['bglive_{}_{}_path'.format(iel, pft_i)],
            temp_val_dict['statv_temp'])
        raster_sum(
            temp_val_dict['statv_temp'], _SV_NODATA,
            temp_val_dict['uptake_below'], _TARGET_NODATA,
            sv_reg['bglive_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)

    # uptake from N fixation into above and belowground live
    if iel == 1:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['uptake_Nfix'], eup_above_iel_path,
                eup_below_iel_path]],
            calc_aboveground_uptake, temp_val_dict['uptake_above'],
            gdal.GDT_Float32, _TARGET_NODATA)
        shutil.copyfile(
            sv_reg['aglive_{}_{}_path'.format(iel, pft_i)],
            temp_val_dict['statv_temp'])
        raster_sum(
            temp_val_dict['statv_temp'], _SV_NODATA,
            temp_val_dict['uptake_above'], _TARGET_NODATA,
            sv_reg['aglive_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['uptake_Nfix'], eup_above_iel_path,
                eup_below_iel_path]],
            calc_belowground_uptake, temp_val_dict['uptake_below'],
            gdal.GDT_Float32, _TARGET_NODATA)
        shutil.copyfile(
            sv_reg['bglive_{}_{}_path'.format(iel, pft_i)],
            temp_val_dict['statv_temp'])
        raster_sum(
            temp_val_dict['statv_temp'], _SV_NODATA,
            temp_val_dict['uptake_below'], _TARGET_NODATA,
            sv_reg['bglive_{}_{}_path'.format(iel, pft_i)], _SV_NODATA)


def calc_nutrient_limitation(return_type):
    """Calculate C, N, and P in new production given nutrient availability."""
    def _nutrlm(
            potenc, rtsh, eavail_1, eavail_2, snfxmx_1, cercrp_max_above_1,
            cercrp_max_below_1, cercrp_max_above_2, cercrp_max_below_2,
            cercrp_min_above_1, cercrp_min_below_1, cercrp_min_above_2,
            cercrp_min_below_2):
        """Calculate new production limited by N and P.

        Growth of new biomass is limited by the availability of N and P.
        Compare nutrient availability to the demand for each nutrient, which
        differs between above- and belowground production. Nutrlm.f

        Parameters:
        potenc (numpy.ndarray): derived, potential production of C calculated
            by root:shoot ratio submodel
        rtsh (numpy.ndarray): derived, root/shoot ratio of new production
        eavail_1 (numpy.ndarray): derived, available N (includes predicted
            symbiotic N fixation)
        eavail_2 (numpy.ndarray): derived, available P
        snfxmx_1 (numpy.ndarray): parameter, maximum symbiotic N fixation rate
        cercrp_max_above_1 (numpy.ndarray): max C/N ratio of new aboveground
            growth
        cercrp_max_below_1 (numpy.ndarray): max C/N ratio of new belowground
            growth
        cercrp_max_above_2 (numpy.ndarray): max C/P ratio of new aboveground
            growth
        cercrp_max_below_2 (numpy.ndarray): max C/P ratio of new belowground
            growth
        cercrp_min_above_1 (numpy.ndarray): min C/N ratio of new aboveground
            growth
        cercrp_min_below_1 (numpy.ndarray): min C/N ratio of new belowground
            growth
        cercrp_min_above_2 (numpy.ndarray): min C/P ratio of new aboveground
            growth
        cercrp_min_below_2 (numpy.ndarray): min C/P ratio of new belowground
            growth

        Returns:
            cprodl, total C production limited by nutrient availability,
                if return_type is 'cprodl'
            eup_above_1, N in new aboveground production, if return_type is
                'eup_above_1'
            eup_below_1, N in new belowground production, if return_type is
                'eup_below_1'
            eup_above_2, P in new aboveground production, if return_type is
                'eup_above_2'
            eup_below_2, P in new belowground production, if return_type is
                'eup_below_2'
            plantNfix, N fixation that actually occurs, if return_type is
                'plantNfix'
        """
        valid_mask = (
            (potenc != _TARGET_NODATA) &
            (rtsh != _TARGET_NODATA) &
            (eavail_1 != _TARGET_NODATA) &
            (eavail_2 != _TARGET_NODATA) &
            (snfxmx_1 != _IC_NODATA) &
            (cercrp_max_above_1 != _TARGET_NODATA) &
            (cercrp_max_below_1 != _TARGET_NODATA) &
            (cercrp_max_above_2 != _TARGET_NODATA) &
            (cercrp_max_below_2 != _TARGET_NODATA) &
            (cercrp_min_above_1 != _TARGET_NODATA) &
            (cercrp_min_below_1 != _TARGET_NODATA) &
            (cercrp_min_above_2 != _TARGET_NODATA) &
            (cercrp_min_below_2 != _TARGET_NODATA))
        cfrac_below = numpy.empty(potenc.shape, dtype=numpy.float32)
        cfrac_below[valid_mask] = (
            rtsh[valid_mask] / (rtsh[valid_mask] + 1.))
        cfrac_above = numpy.empty(potenc.shape, dtype=numpy.float32)
        cfrac_above[valid_mask] = 1. - cfrac_below[valid_mask]

        # maxec is average e/c ratio across aboveground and belowground
        # maxeci is indexed to aboveground only or belowground only
        maxeci_above_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        mineci_above_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        maxeci_below_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        mineci_below_1 = numpy.empty(potenc.shape, dtype=numpy.float32)

        maxeci_above_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        mineci_above_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        maxeci_below_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        mineci_below_2 = numpy.empty(potenc.shape, dtype=numpy.float32)

        maxeci_above_1[valid_mask] = 1. / cercrp_min_above_1[valid_mask]
        mineci_above_1[valid_mask] = 1. / cercrp_max_above_1[valid_mask]
        maxeci_below_1[valid_mask] = 1. / cercrp_min_below_1[valid_mask]
        mineci_below_1[valid_mask] = 1. / cercrp_max_below_1[valid_mask]

        maxeci_above_2[valid_mask] = 1. / cercrp_min_above_2[valid_mask]
        mineci_above_2[valid_mask] = 1. / cercrp_max_above_2[valid_mask]
        maxeci_below_2[valid_mask] = 1. / cercrp_min_below_2[valid_mask]
        mineci_below_2[valid_mask] = 1. / cercrp_max_below_2[valid_mask]

        maxec_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        maxec_1[valid_mask] = (
            cfrac_below[valid_mask] * maxeci_below_1[valid_mask] +
            cfrac_above[valid_mask] * maxeci_above_1[valid_mask])
        maxec_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        maxec_2[valid_mask] = (
            cfrac_below[valid_mask] * maxeci_below_2[valid_mask] +
            cfrac_above[valid_mask] * maxeci_above_2[valid_mask])

        # N/C ratio in new production according to demand and supply
        demand_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        demand_1[valid_mask] = potenc[valid_mask] * maxec_1[valid_mask]

        ecfor_above_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        ecfor_below_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        ecfor_above_1[valid_mask] = (
            mineci_above_1[valid_mask] +
            (maxeci_above_1[valid_mask] - mineci_above_1[valid_mask]) *
            eavail_1[valid_mask] / demand_1[valid_mask])
        ecfor_below_1[valid_mask] = (
            mineci_below_1[valid_mask] +
            (maxeci_below_1[valid_mask] - mineci_below_1[valid_mask]) *
            eavail_1[valid_mask] / demand_1[valid_mask])

        sufficient_mask = ((eavail_1 > demand_1) & valid_mask)
        ecfor_above_1[sufficient_mask] = maxeci_above_1[sufficient_mask]
        ecfor_below_1[sufficient_mask] = maxeci_below_1[sufficient_mask]

        # caculate C production limited by N supply
        c_constrained_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        c_constrained_1[valid_mask] = (
            eavail_1[valid_mask] / (
                cfrac_below[valid_mask] * ecfor_below_1[valid_mask] +
                cfrac_above[valid_mask] * ecfor_above_1[valid_mask]))

        # P/C ratio in new production according to demand and supply
        demand_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        demand_2[valid_mask] = potenc[valid_mask] * maxec_2[valid_mask]

        ecfor_above_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        ecfor_below_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        ecfor_above_2[valid_mask] = (
            mineci_above_2[valid_mask] +
            (maxeci_above_2[valid_mask] - mineci_above_2[valid_mask]) *
            eavail_2[valid_mask] / demand_2[valid_mask])
        ecfor_below_2[valid_mask] = (
            mineci_below_2[valid_mask] +
            (maxeci_below_2[valid_mask] - mineci_below_2[valid_mask]) *
            eavail_2[valid_mask] / demand_2[valid_mask])

        sufficient_mask = ((eavail_2 > demand_2) & valid_mask)
        ecfor_above_2[sufficient_mask] = maxeci_above_2[sufficient_mask]
        ecfor_below_2[sufficient_mask] = maxeci_below_2[sufficient_mask]

        # caculate C production limited by P supply
        c_constrained_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        c_constrained_2[valid_mask] = (
            eavail_2[valid_mask] / (
                cfrac_below[valid_mask] * ecfor_below_2[valid_mask] +
                cfrac_above[valid_mask] * ecfor_above_2[valid_mask]))

        # C production limited by both N and P
        cprodl = numpy.empty(potenc.shape, dtype=numpy.float32)
        cprodl[:] = _TARGET_NODATA
        cprodl[valid_mask] = numpy.minimum(
            c_constrained_1[valid_mask],
            c_constrained_2[valid_mask])
        cprodl[valid_mask] = numpy.minimum(
            cprodl[valid_mask], potenc[valid_mask])

        # N uptake into new production, given limited C production
        eup_above_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        eup_below_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        eup_above_1[:] = _TARGET_NODATA
        eup_below_1[:] = _TARGET_NODATA

        eup_above_1[valid_mask] = (
            cprodl[valid_mask] * cfrac_above[valid_mask] *
            ecfor_above_1[valid_mask])
        eup_below_1[valid_mask] = (
            cprodl[valid_mask] * cfrac_below[valid_mask] *
            ecfor_below_1[valid_mask])

        # P uptake into new production, given limited C production
        eup_above_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        eup_below_2 = numpy.empty(potenc.shape, dtype=numpy.float32)
        eup_above_2[:] = _TARGET_NODATA
        eup_below_2[:] = _TARGET_NODATA

        eup_above_2[valid_mask] = (
            cprodl[valid_mask] * cfrac_above[valid_mask] *
            ecfor_above_2[valid_mask])
        eup_below_2[valid_mask] = (
            cprodl[valid_mask] * cfrac_below[valid_mask] *
            ecfor_below_2[valid_mask])

        # Calculate N fixation that occurs to subsidize needed N supply
        maxNfix = numpy.empty(potenc.shape, dtype=numpy.float32)
        maxNfix[valid_mask] = snfxmx_1[valid_mask] * potenc[valid_mask]

        eprodl_1 = numpy.empty(potenc.shape, dtype=numpy.float32)
        eprodl_1[valid_mask] = (
            eup_above_1[valid_mask] + eup_below_1[valid_mask])
        Nfix_mask = (
            (eprodl_1 - (eavail_1 + maxNfix) > 0.05) &
            valid_mask)
        eprodl_1[Nfix_mask] = eavail_1[Nfix_mask] + maxNfix[Nfix_mask]

        plantNfix = numpy.empty(potenc.shape, dtype=numpy.float32)
        plantNfix[:] = _TARGET_NODATA
        plantNfix[valid_mask] = numpy.maximum(
            eprodl_1[valid_mask] - eavail_1[valid_mask], 0.)

        if return_type == 'cprodl':
            return cprodl
        elif return_type == 'eup_above_1':
            return eup_above_1
        elif return_type == 'eup_below_1':
            return eup_below_1
        elif return_type == 'eup_above_2':
            return eup_above_2
        elif return_type == 'eup_below_2':
            return eup_below_2
        elif return_type == 'plantNfix':
            return plantNfix
    return _nutrlm


def _new_growth(
        pft_id_set, aligned_inputs, site_param_table, veg_trait_table, sv_reg,
        month_reg, current_month):
    """Growth of new aboveground and belowground biomass.

    Add new growth to aboveground and belowground live biomass. C is taken up
    from the atmosphere, while N and P are taken up from the crop storage pool,
    soil mineral N and P content, and symbiotic N fixation. N and P taken up
    from the soil are weighted by the percent cover of the plant functional
    type.

    Parameters:
        pft_id_set (set): set of integers identifying plant functional types
        aligned_inputs (dict): map of key, path pairs indicating paths
            to aligned model inputs, including percent cover of each plant
            functional type
        site_param_table (dict): map of site spatial index to dictionaries
            that contain site-level parameters
        veg_trait_table (dict): map of pft id to dictionaries containing
            plant functional type parameters
        sv_reg (dict): map of key, path pairs giving paths to state variables
            for the current month
        month_reg (dict): map of key, path pairs giving paths to intermediate
            calculated values that are shared between submodels
        current_month (int): month of the year, such that current_month=1
            indicates January

    Modifies:
        The rasters indicated by
            sv_reg['aglivc_<pft>_path'] for each pft
            sv_reg['aglive_1_<pft>_path'] for each pft
            sv_reg['aglive_2_<pft>_path'] for each pft
            sv_reg['bglivc_<pft>_path'] for each pft
            sv_reg['bglive_1_<pft>_path'] for each pft
            sv_reg['bglive_2_<pft>_path'] for each pft
            sv_reg['crpstg_1_<pft>_path'] for each pft
            sv_reg['crpstg_2_<pft>_path'] for each pft
            sv_reg['minerl_{layer}_1_path'] for each soil layer
            sv_reg['minerl_{layer}_2_path'] for each soil layer

    Returns:
        None
    """
    temp_dir = tempfile.mkdtemp(dir=PROCESSING_DIR)
    temp_val_dict = {}
    for val in ['statv_temp']:
        target_path = os.path.join(
            temp_dir, '{}.tif'.format(val))
        temp_val_dict['{}'.format(val)] = target_path
    for val in [
            'availm_1', 'availm_2', 'eavail_1', 'eavail_2', 'potenc',
            'potenc_lim_minerl', 'cprodl', 'eup_above_1', 'eup_below_1',
            'eup_above_2', 'eup_below_2', 'plantNfix']:
        for pft_i in pft_id_set:
            target_path = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
            temp_val_dict['{}_{}'.format(val, pft_i)] = target_path

    param_val_dict = {}
    # site-level parameters
    for val in [
            'favail_1', 'favail_4', 'favail_5', 'favail_6', 'pslsrb',
            'sorpmx']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        param_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
                (site_code, table) in site_param_table.iteritems()])
        pygeoprocessing.reclassify_raster(
            (aligned_inputs['site_index'], 1), site_to_val, target_path,
            gdal.GDT_Float32, _IC_NODATA)
    param_val_dict['favail_2'] = os.path.join(temp_dir, 'favail_2.tif')
    _calc_favail_P(sv_reg, param_val_dict)

    # pft-level parameters
    for pft_i in pft_id_set:
        for val in ['snfxmx_1']:
            target_path = os.path.join(
                temp_dir, '{}_{}.tif'.format(val, pft_i))
            param_val_dict['{}_{}'.format(val, pft_i)] = target_path
            fill_val = veg_trait_table[pft_i][val]
            pygeoprocessing.new_raster_from_base(
                sv_reg['aglivc_{}_path'.format(pft_i)], target_path,
                gdal.GDT_Float32, [_IC_NODATA], fill_value_list=[fill_val])

    for pft_i in pft_id_set:
        if str(current_month) != veg_trait_table[pft_i]['senescence_month']:
            # calculate available nutrients for all pfts prior to
            # performing uptake
            for iel in [1, 2]:
                _calc_avail_mineral_nutrient(
                    veg_trait_table[pft_i], sv_reg, iel,
                    temp_val_dict['availm_{}_{}'.format(iel, pft_i)])
    for pft_i in pft_id_set:
        # growth only occurs in months when senescence not scheduled
        if str(current_month) != veg_trait_table[pft_i]['senescence_month']:
            # calculate available nutrients
            for iel in [1, 2]:
                # eavail_iel, available nutrient
                _calc_available_nutrient(
                    pft_i, iel, veg_trait_table[pft_i], sv_reg,
                    site_param_table, aligned_inputs['site_index'],
                    temp_val_dict['availm_{}_{}'.format(iel, pft_i)],
                    param_val_dict['favail_{}'.format(iel)],
                    month_reg['tgprod_pot_prod_{}'.format(pft_i)],
                    temp_val_dict['eavail_{}_{}'.format(iel, pft_i)])

            # convert from grams of biomass to grams of carbon
            convert_biomass_to_C(
                month_reg['tgprod_pot_prod_{}'.format(pft_i)],
                temp_val_dict['potenc_{}'.format(pft_i)])

            # restrict potential growth by availability of N and P
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['potenc_{}'.format(pft_i)],
                    temp_val_dict['availm_1_{}'.format(pft_i)],
                    temp_val_dict['availm_2_{}'.format(pft_i)],
                    param_val_dict['snfxmx_1_{}'.format(pft_i)]]],
                restrict_potential_growth,
                temp_val_dict['potenc_lim_minerl_{}'.format(pft_i)],
                gdal.GDT_Float32, _TARGET_NODATA)
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['potenc_lim_minerl_{}'.format(pft_i)],
                    month_reg['rtsh_{}'.format(pft_i)],
                    temp_val_dict['eavail_1_{}'.format(pft_i)],
                    temp_val_dict['eavail_2_{}'.format(pft_i)],
                    param_val_dict['snfxmx_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_2_{}'.format(pft_i)]]],
                calc_nutrient_limitation('cprodl'),
                temp_val_dict['cprodl_{}'.format(pft_i)],
                gdal.GDT_Float32, _TARGET_NODATA)
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['potenc_lim_minerl_{}'.format(pft_i)],
                    month_reg['rtsh_{}'.format(pft_i)],
                    temp_val_dict['eavail_1_{}'.format(pft_i)],
                    temp_val_dict['eavail_2_{}'.format(pft_i)],
                    param_val_dict['snfxmx_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_2_{}'.format(pft_i)]]],
                calc_nutrient_limitation('eup_above_1'),
                temp_val_dict['eup_above_1_{}'.format(pft_i)],
                gdal.GDT_Float32, _TARGET_NODATA)
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['potenc_lim_minerl_{}'.format(pft_i)],
                    month_reg['rtsh_{}'.format(pft_i)],
                    temp_val_dict['eavail_1_{}'.format(pft_i)],
                    temp_val_dict['eavail_2_{}'.format(pft_i)],
                    param_val_dict['snfxmx_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_2_{}'.format(pft_i)]]],
                calc_nutrient_limitation('eup_below_1'),
                temp_val_dict['eup_below_1_{}'.format(pft_i)],
                gdal.GDT_Float32, _TARGET_NODATA)
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['potenc_lim_minerl_{}'.format(pft_i)],
                    month_reg['rtsh_{}'.format(pft_i)],
                    temp_val_dict['eavail_1_{}'.format(pft_i)],
                    temp_val_dict['eavail_2_{}'.format(pft_i)],
                    param_val_dict['snfxmx_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_2_{}'.format(pft_i)]]],
                calc_nutrient_limitation('eup_above_2'),
                temp_val_dict['eup_above_2_{}'.format(pft_i)],
                gdal.GDT_Float32, _TARGET_NODATA)
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['potenc_lim_minerl_{}'.format(pft_i)],
                    month_reg['rtsh_{}'.format(pft_i)],
                    temp_val_dict['eavail_1_{}'.format(pft_i)],
                    temp_val_dict['eavail_2_{}'.format(pft_i)],
                    param_val_dict['snfxmx_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_2_{}'.format(pft_i)]]],
                calc_nutrient_limitation('eup_below_2'),
                temp_val_dict['eup_below_2_{}'.format(pft_i)],
                gdal.GDT_Float32, _TARGET_NODATA)
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['potenc_lim_minerl_{}'.format(pft_i)],
                    month_reg['rtsh_{}'.format(pft_i)],
                    temp_val_dict['eavail_1_{}'.format(pft_i)],
                    temp_val_dict['eavail_2_{}'.format(pft_i)],
                    param_val_dict['snfxmx_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_max_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_max_below_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_1_{}'.format(pft_i)],
                    month_reg['cercrp_min_above_2_{}'.format(pft_i)],
                    month_reg['cercrp_min_below_2_{}'.format(pft_i)]]],
                calc_nutrient_limitation('plantNfix'),
                temp_val_dict['plantNfix_{}'.format(pft_i)],
                gdal.GDT_Float32, _TARGET_NODATA)

            # uptake of C into new aboveground production
            shutil.copyfile(
                sv_reg['aglivc_{}_path'.format(pft_i)],
                temp_val_dict['statv_temp'])
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['statv_temp'],
                    temp_val_dict['cprodl_{}'.format(pft_i)],
                    month_reg['rtsh_{}'.format(pft_i)]]],
                c_uptake_aboveground, sv_reg['aglivc_{}_path'.format(pft_i)],
                gdal.GDT_Float32, _SV_NODATA)

            # uptake of C into new belowground production
            shutil.copyfile(
                sv_reg['bglivc_{}_path'.format(pft_i)],
                temp_val_dict['statv_temp'])
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [
                    temp_val_dict['statv_temp'],
                    temp_val_dict['cprodl_{}'.format(pft_i)],
                    month_reg['rtsh_{}'.format(pft_i)]]],
                c_uptake_belowground, sv_reg['bglivc_{}_path'.format(pft_i)],
                gdal.GDT_Float32, _SV_NODATA)

            # uptake of N and P into new above- and belowground production
            for iel in [1, 2]:
                nutrient_uptake(
                    iel, int(veg_trait_table[pft_i]['nlaypg']),
                    aligned_inputs['pft_{}'.format(pft_i)],
                    temp_val_dict['eup_above_{}_{}'.format(iel, pft_i)],
                    temp_val_dict['eup_below_{}_{}'.format(iel, pft_i)],
                    temp_val_dict['plantNfix_{}'.format(pft_i)],
                    temp_val_dict['availm_{}_{}'.format(iel, pft_i)],
                    temp_val_dict['eavail_{}_{}'.format(iel, pft_i)],
                    sv_reg, pft_i, param_val_dict['pslsrb'],
                    param_val_dict['sorpmx'])