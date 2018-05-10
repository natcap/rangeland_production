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
import numpy as np

import pygeoprocessing
from natcap.invest import utils

LOGGER = logging.getLogger('natcap.invest.forage')

# we only have these types of soils
SOIL_TYPE_LIST = ['clay', 'silt', 'sand']
# state variables that are updated at each monthly timestep
# state variables property of the site
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
# state variables property of PFT
_PFT_STATE_VARIABLES = ['aglivc', 'bglivc', 'stdedc', 'aglive_1', 'bglive_1',
    'stdede_1', 'aglive_2', 'bglive_2', 'stdede_2']
# intermediate parameters that do not change between timesteps
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

# Target nodata is for general rasters that are positive, and _IC_NODATA are
# for rasters that are any range
_TARGET_NODATA = -1.0
_IC_NODATA = np.finfo('float32').min

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
        args['site_param_table'] string: path to csv file giving site
            parameters. This file must contain a column named "site" that
            contains unique integers. These integer values correspond to site
            type identifiers which are values in the site parameter spatial
            index raster. Other required fields for this table are site and
            "fixed" parameters from the Century model, i.e., the parameters
            in the Century input files site.100 and fix.100.
        args['site_param_spatial_index_path'] string: path to a raster file
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
        args['animal_trait_path'] string: path to csv file giving animal traits
            for each animal type - number - duration combination. This table
            must contain a column named "animal_id" that contains unique
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
    for month_index in xrange(int(args['n_months'])):
        month_i = (starting_month + month_index - 1) % 12 + 1
        temperature_month_set.add(month_i)
        year = starting_year +  (starting_month + month_index - 1) // 12
        precip_path = args['monthly_precip_path_pattern'].replace(
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
    # atmospheric deposition of N and S from annual precipitation
    n_precip_months = int(args['n_months'])
    if n_precip_months < 12:
        month_index = int(args['n_months'])
        while month_index <= 12:
            month_i = (starting_month + month_index - 1) % 12 + 1
            year = starting_year +  (starting_month + month_index - 1) // 12
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
            base_align_raster_path_id_map['%s_temp_%d' \
                % (substring, month_i)] = monthly_temp_path
            if not os.path.exists(monthly_temp_path):
                missing_min_temperature_path_list.append(monthly_temp_path)
        if missing_temperature_path_list:
            raise ValueError(
            "Couldn't find the following temperature raster paths given the " +
            "pattern: %s\n\t" % args['%s_temp_path_pattern' % substring] +
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
    base_align_raster_path_id_map['site_index'] = \
        args['site_param_spatial_index_path']
    assert (pygeoprocessing.get_raster_info(
        args['site_param_spatial_index_path'])['n_bands'] == 1), """Site
        spatial index raster must contain only one band"""
    assert (pygeoprocessing.get_raster_info(
        args['site_param_spatial_index_path'])['datatype'] in
        [1, 2, 3, 4, 5]), "Site spatial index raster must be integer type"
    # get unique values in site param raster
    site_index_set = set()
    for offset_map, raster_block in pygeoprocessing.iterblocks(
            args['site_param_spatial_index_path']):
        site_index_set.update(np.unique(raster_block))
    site_nodata = pygeoprocessing.get_raster_info(
        args['site_param_spatial_index_path'])['nodata'][0]
    if site_nodata in site_index_set:
        site_index_set.remove(site_nodata)
    site_param_table = utils.build_lookup_from_csv(args['site_param_table'],
                                                  'site')
    missing_site_index_list = list(site_index_set.difference(
        set(site_param_table.keys())))
    if missing_site_index_list:
        raise ValueError(
            "Couldn't find parameter values for the following site " +
            "indices: %s\n\t" + ", ".join(missing_site_index_list))
            
    # make sure veg traits exist for each pft raster
    pft_dir = os.path.dirname(args['veg_spatial_composition_path_pattern'])
    pft_basename = os.path.basename(
        args['veg_spatial_composition_path_pattern'])
    files = [f for f in os.listdir(pft_dir) if os.path.isfile(
             os.path.join(pft_dir, f))]
    pft_regex = re.compile(pft_basename.replace('<PFT>', '(\d+)'))
    pft_matches = [
        m for m in [pft_regex.search(f) for f in files] if m is not None]
    pft_id_set = set([int(m.group(1)) for m in pft_matches])
    for pft_i in pft_id_set:
        pft_path = args['veg_spatial_composition_path_pattern'].replace(
            '<PFT>', '%d' % pft_i)
        base_align_raster_path_id_map['pft_%d' % pft_i] = pft_path
    veg_trait_table = utils.build_lookup_from_csv(args['veg_trait_path'],
                                                  'PFT')
    missing_pft_trait_list = pft_id_set.difference(
        set(veg_trait_table.keys()))
    if missing_pft_trait_list:
        raise ValueError(
            "Couldn't find trait values for the following plant functional " +
            "types: %s\n\t" + ", ".join(missing_pft_trait_list))
    
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
    
    anim_trait_table = utils.build_lookup_from_csv(args['animal_trait_path'],
                                                   'animal_id')
    missing_animal_trait_list = set(anim_id_list).difference(
        set(anim_trait_table.keys()))
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
    source_input_path_list = [base_align_raster_path_id_map[k] for k in
        sorted(base_align_raster_path_id_map.iterkeys())]
    aligned_input_path_list = [aligned_inputs[k] for k in
        sorted(aligned_inputs.iterkeys())]
    
    file_suffix = utils.make_suffix_string(args, 'results_suffix')
    
    # if initial conditions are supplied, use them to initialize all state
    # variables
    if args['initial_conditions_dir']:
        LOGGER.info("setting initial conditions from this directory: %s",
                    args['initial_conditions_dir'])
        
        # check that all necessary state variables are supplied
        missing_initial_values = []
        # align initial state variables to resampled inputs
        resample_initial_path_map = {}
        for sv in _SITE_STATE_VARIABLE_FILES.keys():
            sv_path = os.path.join(args['initial_conditions_dir'],
                _SITE_STATE_VARIABLE_FILES[sv])
            resample_initial_path_map[sv] = sv_path
            if not os.path.exists(sv_path):
                missing_initial_values.append(sv_path)
        for pft_index in pft_id_set:
            for sv in _PFT_STATE_VARIABLES:
                sv_key = '{}_{}_path'.format(sv, pft_index)
                sv_path = os.path.join(args['initial_conditions_dir'],
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
        aligned_initial_path_map = dict([(key, os.path.join(sv_dir,
            os.path.basename(path))) for key, path in
            resample_initial_path_map.iteritems()])
        initial_path_list = ([resample_initial_path_map[k] for k in
            sorted(resample_initial_path_map.iterkeys())] +
            source_input_path_list)
        aligned_initial_path_list = ([aligned_initial_path_map[k] for k in
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
    
    ## Initialization
    # calculate persistent intermediate parameters that do not change during
    # the simulation
    persist_param_dir = os.path.join(args['workspace_dir'], 
        'intermediate_parameters')
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
    _structural_ratios(aligned_inputs['site_index'], site_param_table,
        sv_reg, pp_reg)
    
    ## Main simulation loop
    # for each step in the simulation
    for month_index in xrange(int(args['n_months'])):
        # make new folders for state variables during this step
        sv_dir = os.path.join(args['workspace_dir'],
            'state_variables_m%d' % month_index)
        utils.make_directories([sv_dir])
        
        # track state variables from previous step
        prev_sv_reg = sv_reg
        sv_reg = utils.build_file_registry(
            [(_SITE_STATE_VARIABLE_FILES, sv_dir),
            (pft_sv_dict, sv_dir)], file_suffix)
        
        # update state variables from previous month
        LOGGER.info("Main simulation loop: month %d of %d" % (month_index,
            int(args['n_months'])))
            
def _afiel_awilt(site_index_path, site_param_table, som1c_2_path,
                 som2c_2_path, som3c_path, sand_path, silt_path, clay_path,
                 bulk_d_path, pp_reg):
    """Calculate field capacity and wilting point for each soil layer.
    Computations based on Gupta and Larson 1979, 'Estimating soil and water
    retention characteristics from particle size distribution, organic 
    matter percent and bulk density'. Water Resources Research 15:1633.
    Field capacity is calculated for -0.33 bar; wilting point is
    calculated for water content at -15 bars."""
    
    # temporary intermediate rasters for this calculation
    temp_dir = tempfile.mkdtemp()
    edepth_path = os.path.join(temp_dir, 'edepth.tif')
    ompc_path = os.path.join(temp_dir, 'ompc.tif')
    
    def calc_ompc(som1c_2, som2c_2, som3c, bulkd, edepth):
        """Calculate organic matter as the sum of soil organic matter,
        weighted by bulk density. From line 222, Prelim.f"""
        ompc = np.empty(som1c_2.shape)
        ompc[:] = _TARGET_NODATA
        valid_mask = ((som1c_2 > 0) & (som2c_2 > 0) & (som3c > 0) &
            (bulkd > 0))
        ompc[valid_mask] = (som1c_2[valid_mask] + som2c_2[valid_mask] +
            som3c[valid_mask]) * 1.724 / (10000. * bulkd[valid_mask] *
            edepth[valid_mask])
        return ompc

    def calc_afiel(sand, silt, clay, ompc, bulkd):
        """Calculate field capacity."""
        afiel = np.empty(sand.shape)
        afiel[:] = _TARGET_NODATA
        valid_mask = ompc != _TARGET_NODATA
        afiel[valid_mask] = (0.3075 * sand[valid_mask] +
            0.5886 * silt[valid_mask] + 0.8039 * clay[valid_mask] + 
            2.208E-03 * ompc[valid_mask] + -0.1434 * bulkd[valid_mask])
        return afiel
    
    def calc_awilt(sand, silt, clay, ompc, bulkd):
        """Calculate wilting point."""
        awilt = np.empty(sand.shape)
        awilt[:] = _TARGET_NODATA
        valid_mask = ompc != _TARGET_NODATA
        awilt[valid_mask] = (-0.0059 * sand[valid_mask]
            + 0.1142 * silt[valid_mask] + 0.5766 * clay[valid_mask] +
            2.228E-03 * ompc[valid_mask] + 0.02671 * bulkd[valid_mask])
        return awilt
    
    def decrement_ompc(ompc):
        return ompc * 0.85
    
    # temporary raster from site parameter 'edepth'
    site_to_edepth = dict(
        [(site_code, float(table['edepth'])) for
         (site_code, table) in site_param_table.items()])

    pygeoprocessing.reclassify_raster(
        (site_index_path, 1), site_to_edepth, edepth_path, gdal.GDT_Float32,
        _TARGET_NODATA)
    
    # estimate total soil organic matter
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [som1c_2_path, som2c_2_path, som3c_path,
            bulk_d_path, edepth_path]],
        calc_ompc, ompc_path, gdal.GDT_Float32, _TARGET_NODATA)
    
    # calculate field capacity and wilting point for each soil layer,
    # decreasing organic matter content by 85% with each layer
    for lyr in xrange(1, 10):
        afiel_path = pp_reg['afiel_{}_path'.format(lyr)]
        awilt_path = pp_reg['awilt_{}_path'.format(lyr)]
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [sand_path, silt_path, clay_path,
                ompc_path, bulk_d_path]],
            calc_afiel, afiel_path, gdal.GDT_Float32, _TARGET_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [sand_path, silt_path, clay_path,
                ompc_path, bulk_d_path]],
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
    """Calculate persistent parameters that do not change over the course of
    the simulation."""
    
    # temporary intermediate rasters for these calculations
    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {}
    for val in['peftxa', 'peftxb', 'p1co2a_2', 'p1co2b_2', 'ps1s3_1',
               'ps1s3_2', 'ps2s3_1', 'ps2s3_2', 'omlech_1', 'omlech_2']:
        target_path = os.path.join(temp_dir, '{}.tif'.format(val))
        temp_val_dict[val] = target_path
        site_to_val = dict(
            [(site_code, float(table[val])) for
            (site_code, table) in site_param_table.items()])
        pygeoprocessing.reclassify_raster(
            (site_index_path, 1), site_to_val, target_path, gdal.GDT_Float32,
            _TARGET_NODATA)
    
    def calc_eftext(peftxa, peftxb, sand):
        """Calculate the effect of soil texture on soil microbe decomposition
        rate. Line 359 Prelim.f"""
        eftext = np.empty(sand.shape)
        eftext[:] = _IC_NODATA
        valid_mask = sand > 0
        eftext[valid_mask] = peftxa[valid_mask] + (
            peftxb[valid_mask] * sand[valid_mask])
        return eftext
    
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [temp_val_dict['peftxa'],
            temp_val_dict['peftxb'], sand_path]],
        calc_eftext, pp_reg['eftext_path'], gdal.GDT_Float32, _IC_NODATA)
    
    def calc_p1co2_2(p1co2a_2, p1co2b_2, sand):
        """Calculate the fraction of C lost to CO2 during decomposition from
        som1c_2. Line 366 Prelim.f"""
        p1co2_2 = np.empty(sand.shape)
        p1co2_2[:] = _IC_NODATA
        valid_mask = sand > 0
        p1co2_2[valid_mask] = p1co2a_2[valid_mask] + (
            p1co2b_2[valid_mask] * sand[valid_mask])
        return p1co2_2
    
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [temp_val_dict['p1co2a_2'],
            temp_val_dict['p1co2b_2'], sand_path]],
        calc_p1co2_2, pp_reg['p1co2_2_path'], gdal.GDT_Float32, _IC_NODATA)
        
    def calc_fps1s3(ps1s3_1, ps1s3_2, clay):
        """Calculate the effect of clay content on decomposition from
        som1c_2.  Line 370 Prelim.f"""
        fps1s3 = np.empty(clay.shape)
        fps1s3[:] = _IC_NODATA
        valid_mask = clay > 0
        fps1s3[valid_mask] = ps1s3_1[valid_mask] + (
            ps1s3_2[valid_mask] * clay[valid_mask])
        return fps1s3
    
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [temp_val_dict['ps1s3_1'],
            temp_val_dict['ps1s3_2'], clay_path]],
        calc_fps1s3, pp_reg['fps1s3_path'], gdal.GDT_Float32, _IC_NODATA)
            
    def calc_fps2s3(ps2s3_1, ps2s3_2, clay):
        """Calculate the effect of clay content on decomposition from
        som2c_2.  Line 371 Prelim.f"""
        fps2s3 = np.empty(clay.shape)
        fps2s3[:] = _IC_NODATA
        valid_mask = clay > 0
        fps2s3[valid_mask] = ps2s3_1[valid_mask] + (
            ps2s3_2[valid_mask] * clay[valid_mask])
        return fps2s3
    
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [temp_val_dict['ps2s3_1'],
            temp_val_dict['ps2s3_2'], clay_path]],
        calc_fps2s3, pp_reg['fps2s3_path'], gdal.GDT_Float32, _IC_NODATA)
    
    def calc_orglch(omlech_1, omlech_2, sand):
        """Calculate the effect of sand content on organic leaching from
        soil.  Line 110 Predec.f"""
        orglch = np.empty(sand.shape)
        orglch[:] = _IC_NODATA
        valid_mask = sand > 0
        orglch[valid_mask] = omlech_1[valid_mask] + (
            omlech_2[valid_mask] * sand[valid_mask])
        return orglch
    
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [temp_val_dict['omlech_1'],
            temp_val_dict['omlech_2'], sand_path]],
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
    
    # temporary intermediate rasters for these calculations
    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {}
    for iel in [1, 2]:
        for val in['pcemic1_2', 'pcemic1_1', 'pcemic1_3', 'pcemic2_2',
                   'pcemic2_1', 'pcemic2_3', 'rad1p_1', 'rad1p_2',
                   'rad1p_3', 'varat1_1', 'varat22_1']:
            target_path = os.path.join(temp_dir, '{}_{}.tif'.format(val, iel))
            temp_val_dict['{}_{}'.format(val, iel)] = target_path
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
                temp_val_dict['pcemic1_2_{}'.format(iel)],
                temp_val_dict['pcemic1_1_{}'.format(iel)],
                temp_val_dict['pcemic1_3_{}'.format(iel)],
                sv_reg['struce_1_{}_path'.format(iel)],
                sv_reg['strucc_1_path']]],
            calc_rnewas_som1, pp_reg['rnewas_1_{}_path'.format(iel)],
            gdal.GDT_Float32, _TARGET_NODATA)
        # calculate rnewas_iel_2 - aboveground material to SOM2
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                temp_val_dict['pcemic2_2_{}'.format(iel)],
                temp_val_dict['pcemic2_1_{}'.format(iel)],
                temp_val_dict['pcemic2_3_{}'.format(iel)],
                sv_reg['struce_1_{}_path'.format(iel)],
                sv_reg['strucc_1_path'],
                temp_val_dict['rad1p_1_{}'.format(iel)],
                temp_val_dict['rad1p_2_{}'.format(iel)],
                temp_val_dict['rad1p_3_{}'.format(iel)],
                temp_val_dict['pcemic1_2_{}'.format(iel)],
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
    