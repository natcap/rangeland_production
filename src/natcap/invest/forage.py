"""InVEST Forage module.

InVEST Forage model developed from this design doc:
https://docs.google.com/document/d/10oJo43buEdJkFTZ0wYaW00EagSzs1oM7g_lBUc8URMI/edit#
"""
import os
import logging

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
    }

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
            year. The model expects to find a precipitation file input for
            every month of the simulation, so there should be as many
            precipitation input files as `n_months`. The <month> value in
            input files must be two digits.
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
    base_align_raster_path_id_map['bulk_d'] = args['bulk_density_path']
    base_align_raster_path_id_map['ph'] = args['ph_path']

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
    aligned_raster_path_id_map = dict([(key, os.path.join(
        aligned_raster_dir, 'aligned_%s' % os.path.basename(path)))
        for key, path in base_align_raster_path_id_map.iteritems()])

    # align all the base inputs to be the minimum known pixel size and to
    # only extend over their combined intersections
    source_path_list = [base_align_raster_path_id_map[k] for k in
        sorted(base_align_raster_path_id_map.iterkeys())]
    aligned_path_list = [aligned_raster_path_id_map[k] for k in
        sorted(aligned_raster_path_id_map.iterkeys())]
    LOGGER.info("aligning base raster inputs")
    pygeoprocessing.align_and_resize_raster_stack(
        source_path_list,aligned_path_list,
        ['nearest'] * len(aligned_raster_path_id_map),
        target_pixel_size, 'intersection')
    
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
                    
        # align and resample initialization rasters to match aligned inputs
        sv_dir = os.path.join(args['workspace_dir'], 'state_variables_m-1')
        template_aligned_raster_path = [k for k in 
            aligned_raster_path_id_map.itervalues()][0]
        aligned_initial_path_map = dict([(key, os.path.join(sv_dir,
            os.path.basename(path))) for key, path in
            resample_initial_path_map.iteritems()])
        initial_path_list = [resample_initial_path_map[k] for k in
            sorted(resample_initial_path_map.iterkeys())]
        aligned_initial_path_list = [aligned_initial_path_map[k] for k in
            sorted(aligned_initial_path_map.iterkeys())]
        # insert an aligned input to use as bounding box
        initial_path_list.insert(0, template_aligned_raster_path)
        aligned_initial_path_list.insert(0, os.path.join(sv_dir,
            'align_template.tif'))
        pygeoprocessing.align_and_resize_raster_stack(
            initial_path_list, aligned_initial_path_list,
            ['nearest'] * len(initial_path_list),
            target_pixel_size, 'intersection', raster_align_index=0)
        sv_reg = aligned_initial_path_map
    else:
        LOGGER.info("initial conditions not supplied")
        # TODO add spin-up or initialization with Burke's equations
        raise ValueError("Initial conditions must be supplied")
    
    ## Initialization
    # calculate persistent intermediate parameters that do not change during
    # the simulation
    # make folder for these persistent intermediate parameters
    persist_param_dir = os.path.join(args['workspace_dir'], 
        'intermediate_parameters')
    utils.make_directories([persist_param_dir])
    pp_reg = utils.build_file_registry(
            [(_PERSISTENT_PARAMS_FILES, persist_param_dir)], file_suffix)
    # TODO calculate persistent params
    
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