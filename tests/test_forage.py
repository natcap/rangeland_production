"""InVEST forage model tests."""

import unittest
import tempfile
import shutil
import os
import sys
import math

import numpy
from osgeo import osr
from osgeo import gdal

import pygeoprocessing

SAMPLE_DATA = "C:/Users/ginge/Documents/NatCap/sample_inputs"
REGRESSION_DATA = "C:/Users/ginge/Documents/NatCap/regression_test_data"
PROCESSING_DIR = None

_TARGET_NODATA = -1.0
_IC_NODATA = numpy.finfo('float32').min
_SV_NODATA = -1.0

NROWS = 3
NCOLS = 3

numpy.random.seed(100)


def create_random_raster(
        target_path, lower_bound, upper_bound, nrows=NROWS, ncols=NCOLS):
    """Create a small raster of random floats.

    The raster will have nrows rows and ncols columns and will be in the
    unprojected coordinate system WGS 1984. The values in the raster
    will be between `lower_bound` (included) and `upper_bound`
    (excluded).

    Parameters:
        target_path (string): path to result raster
        lower_bound (float): lower limit of range of random values
            (included)
        upper_bound (float): upper limit of range of random values
            (excluded)

    Returns:
        None
    """
    geotransform = [0, 0.0001, 0, 44.5, 0, 0.0001]
    n_bands = 1
    datatype = gdal.GDT_Float32
    projection = osr.SpatialReference()
    projection.SetWellKnownGeogCS('WGS84')
    driver = gdal.GetDriverByName('GTiff')
    target_raster = driver.Create(
        target_path.encode('utf-8'), ncols, nrows, n_bands,
        datatype)
    target_raster.SetProjection(projection.ExportToWkt())
    target_raster.SetGeoTransform(geotransform)
    target_band = target_raster.GetRasterBand(1)
    target_band.SetNoDataValue(_TARGET_NODATA)

    random_array = numpy.random.uniform(
        lower_bound, upper_bound, (nrows, ncols))
    target_band.WriteArray(random_array)
    target_raster = None


def create_complementary_raster(
        raster1_path, raster2_path, result_raster_path):
    """Create a raster to sum inputs to 1.

    The sum of the two input rasters and the result raster will be
    1.

    Parameters:
        raster1_path (string): path to raster of floast between 0 and 1
        raster2_path (string): path to raster of floast between 0 and 1
        result_raster_path (string): path to result raster

    Modifies:
        the raster indicated by `result_raster_path

    Returns:
        None
    """
    def complement_op(input1, input2):
        """Generate an array that adds to 1 with input1 and input2."""
        result_raster = numpy.empty(input1.shape, dtype=numpy.float32)
        result_raster[:] = _TARGET_NODATA
        valid_mask = (
            (input1 != _TARGET_NODATA)
            & (input2 != _TARGET_NODATA))

        result_raster[valid_mask] = (
            1. - (input1[valid_mask] + input2[valid_mask]))
        return result_raster

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [raster1_path, raster2_path]],
        complement_op, result_raster_path, gdal.GDT_Float32,
        _TARGET_NODATA)


def insert_nodata_values_into_raster(target_raster, nodata_value):
    """Insert nodata at arbitrary locations in `target_raster`."""
    def insert_op(prior_copy):
        modified_copy = prior_copy
        if (prior_copy.shape[0] * prior_copy.shape[1]) == 1:
            n_vals = 1
        else:
            n_vals = numpy.random.randint(
                1, (prior_copy.shape[0] * prior_copy.shape[1]))
        insertions = 0
        while insertions < n_vals:
            row = numpy.random.randint(0, prior_copy.shape[0])
            col = numpy.random.randint(0, prior_copy.shape[1])
            modified_copy[row, col] = nodata_value
            insertions += 1
        return modified_copy

    prior_copy = os.path.join(
        os.path.dirname(target_raster), 'prior_to_insert_nodata.tif')
    shutil.copyfile(target_raster, prior_copy)

    pygeoprocessing.raster_calculator(
        [(prior_copy, 1)],
        insert_op, target_raster,
        gdal.GDT_Float32, nodata_value)

    os.remove(prior_copy)


def create_constant_raster(target_path, fill_value):
    """Create a single-pixel raster with value `fill_value`."""
    geotransform = [0, 1, 0, 44.5, 0, 1]
    n_cols = 1
    n_rows = 1
    n_bands = 1
    datatype = gdal.GDT_Float32
    projection = osr.SpatialReference()
    projection.SetWellKnownGeogCS('WGS84')
    driver = gdal.GetDriverByName('GTiff')
    target_raster = driver.Create(
        target_path.encode('utf-8'), n_cols, n_rows, n_bands,
        datatype)
    target_raster.SetProjection(projection.ExportToWkt())
    target_raster.SetGeoTransform(geotransform)
    target_band = target_raster.GetRasterBand(1)
    target_band.SetNoDataValue(_TARGET_NODATA)
    target_band.Fill(fill_value)
    target_raster = None


def insert_nodata_values_into_array(target_array, nodata_value):
    """Insert nodata at arbitrary locations in `target_array`."""
    modified_array = target_array
    n_vals = numpy.random.randint(
        0, (target_array.shape[0] * target_array.shape[1]))
    insertions = 0
    while insertions < n_vals:
        row = numpy.random.randint(0, target_array.shape[0])
        col = numpy.random.randint(0, target_array.shape[1])
        modified_array[row, col] = nodata_value
        insertions += 1
    return modified_array


def monthly_N_fixation_point(
        precip, annual_precip, baseNdep, epnfs_2, prev_minerl_1_1):
    """Add monthly N fixation to surface mineral N pool.

    Monthly N fixation is calculated from annual N deposition according to
    the ratio of monthly precipitation to annual precipitation.

    Parameters:
        precip (float): input, monthly precipitation
        annual_precip (float): derived, annual precipitation
        baseNdep (float): derived, annual atmospheric N deposition
        epnfs_2 (float): parameter, intercept of regression
            predicting N deposition from annual precipitation
        prev_minerl_1_1 (float): state variable, mineral N in the
            surface layer in previous month

    Returns:
        minerl_1_1, updated mineral N in the surface layer
    """
    wdfxm = (
        baseNdep * (precip / annual_precip) + epnfs_2 *
        min(annual_precip, 100.) * (precip / annual_precip))
    minerl_1_1 = prev_minerl_1_1 + wdfxm
    return minerl_1_1


def rprpet_point(pet, snowmelt, avh2o_3, precip):
    """Calculate the ratio of precipitation to ref evapotranspiration.

    The ratio of precipitation or snowmelt to reference
    evapotranspiration influences agdefac and bgdefac, the above- and
    belowground decomposition factors.

    Parameters:
        pet (float): derived, reference evapotranspiration
        snowmelt (float): derived, snowmelt occuring this month
        avh2o_3 (float): derived, moisture in top two soil layers
        precip (float): input, precipitation for this month

    Returns:
        rprpet, the ratio of precipitation or snowmelt to reference
            evapotranspiration
    """
    if snowmelt > 0:
        rprpet = snowmelt / pet
    else:
        rprpet = (avh2o_3 + precip) / pet
    return rprpet


def defac_point(
        snow, min_temp, max_temp, rprpet, teff_1, teff_2, teff_3, teff_4):
    """Point-based version of `calc_defac`.

    The decomposition factor reflects the influence of soil temperature and
    moisture on decomposition. Lines 151-200, Cycle.f.

    Parameters:
        snow (float): standing snowpack
        min_temp (float): average minimum temperature for the month
        max_temp (float): average maximum temperature for the month
        rprpet (float): ratio of precipitation or snowmelt to
            reference evapotranspiration
        teff_1 (float): x location of inflection point for
            calculating the effect of soil temperature on decomposition
            factor
        teff_2 (float): y location of inflection point for
            calculating the effect of soil temperature on decomposition
            factor
        teff_3 (float): step size for calculating the effect
            of soil temperature on decomposition factor
        teff_4 (float): lope of the line at the inflection
            point, for calculating the effect of soil temperature on
            decomposition factor

    Returns:
        defac, aboveground and belowground decomposition factor
    """
    if rprpet > 9:
        agwfunc = 1
    else:
        agwfunc = 1. / (1 + 30 * math.exp(-8.5 * rprpet))
    if snow > 0:
        stemp = 0
    else:
        stemp = (min_temp + max_temp) / 2.
    tfunc = max(
        0.01, (teff_2 + (teff_3 / math.pi) * numpy.arctan(math.pi *
            teff_4 * (stemp - teff_1))) /
        (teff_2 + (teff_3 / math.pi) * numpy.arctan(math.pi *
            teff_4 * (30.0 - teff_1))))
    defac = max(0, tfunc * agwfunc)
    return defac


def calc_anerb_point(
        rprpet, pevap, drain, aneref_1, aneref_2, aneref_3):
    """Calculate effect of soil anaerobic conditions on decomposition.

    The impact of soil anaerobic conditions on decomposition is
    calculated from soil moisture and reference evapotranspiration.
    Anerob.f.

    Parameters:
        rprpet (float): ratio of precipitation or snowmelt to
            reference evapotranspiration
        pevap (float): reference evapotranspiration
        drain (float): the fraction of excess water lost by
            drainage. Indicates whether a soil is sensitive for
            anaerobiosis (drain = 0) or not (drain = 1)
        aneref_1 (float): value of rprpet below which there
            is no negative impact of soil anaerobic conditions on
            decomposition
        aneref_2 (float): value of rprpet above which there
            is maximum negative impact of soil anaerobic conditions on
            decomposition
        aneref_3 (float): minimum value of the impact of
            soil anaerobic conditions on decomposition

    Returns:
        anerb, the effect of soil anaerobic conditions on decomposition
    """
    anerb = 1
    if rprpet > aneref_1:
        xh2o = (rprpet - aneref_1) * pevap * (1. - drain)
        if xh2o > 0:
            newrat = aneref_1 + (xh2o / pevap)
            slope = (1. - aneref_3) / (aneref_1 - aneref_2)
            anerb = 1. + slope * (newrat - aneref_1)
        anerb = max(anerb, aneref_3)
    return anerb


def bgdrat_point(aminrl, varat_1_iel, varat_2_iel, varat_3_iel):
    """Calculate required C/iel ratio for belowground decomposition.

    When belowground material decomposes, its nutrient content is
    compared to this ratio to check whether nutrient content is
    sufficiently high to allow decomposition. This ratio is calculated at
    each decomposition time step.

    Parameters:
        aminrl (float): mineral <iel> (N or P) in top soil layer, averaged
            across decomposition time steps
        varat_1_iel (float): parameter, maximum C/iel ratio
        varat_2_iel (float): parameter, minimum C/iel ratio
        varat_3_iel (float): parameter, amount of iel present when minimum
            ratio applies

    Returns:
        bgdrat, the required C/iel ratio for decomposition
    """
    if aminrl <= 0:
        bgdrat = varat_1_iel
    elif aminrl > varat_3_iel:
        bgdrat = varat_2_iel
    else:
        bgdrat = (
            (1. - aminrl / varat_3_iel) * (varat_1_iel - varat_2_iel) +
            varat_2_iel)
    return bgdrat


def esched_point(return_type):
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
        """Calculate the flow of one element to accompany decomp of C.

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
            cflow: total C that is decomposing from box A to box B
            tca: C in donating stock, i.e. box A
            rcetob: required ratio of C/iel in the receiving stock
            anps: iel (N or P) in the donating stock
            labile: mineral iel (N or P)

        Returns:
            material_leaving_a, the amount of material leaving box A, if
                return_type is 'material_leaving_a'
            material_arriving_b, the amount of material arriving in box B,
                if return_type is 'material_arriving_b'
            mnrflo, flow in or out of mineral pool, if return_type is
                'mineral_flow'
        """
        outofa = anps * (cflow/tca)
        if (cflow/outofa > rcetob):
            # immobilization occurs
            immflo = cflow/rcetob - outofa
            if ((labile - immflo) > 0):
                material_leaving_a = outofa  # outofa flows from anps to bnps
                material_arriving_b = outofa + immflo
                mnrflo = -immflo  # immflo flows from mineral to bnps
            else:
                mnrflo = 0  # no flow from box A to B, nothing moves
                material_leaving_a = 0
                material_arriving_b = 0
        else:
            # mineralization
            atob = cflow/rcetob
            material_leaving_a = outofa
            material_arriving_b = atob  # atob flows from anps to bnps
            mnrflo = outofa - atob  # the rest of material leaving box A
                                    # goes to mineral
        if return_type == 'material_leaving_a':
            return material_leaving_a
        elif return_type == 'material_arriving_b':
            return material_arriving_b
        elif return_type == 'mineral_flow':
            return mnrflo
    return _esched


def declig_point(return_type):
    """Point implementation of decomposition of structural material.

    Track the decomposition of structural material (i.e., material containing
    lignin) into SOM2 and SOM1.

    Returns:
        the function `_declig`
    """
    def _declig(
            aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr, tcflow,
            struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
            rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2):
        """Decomposition of material containing lignin into SOM2 and SOM1.

        This function is called when surface structural C (STRUCC_1)
        decomposes, and again when soil structural C (STRUCC_2) decomposes.
        Declig.f

        Parameters:
            aminrl_1: mineral N averaged across the 4-times-per-monthly
                decomposition timesteps
            aminrl_2: mineral P averaged across the 4-times-per-monthly
                decomposition timesteps
            ligcon: lignin content of the decomposing material
            rsplig: co2 loss with decomposition to SOM2
            ps1co2_lyr: co2 loss with decomposition to SOM1
            strucc_lyr: structural C in lyr that is decomposing
            tcflow: the total amount of C flowing out of strucc_lyr
            struce_lyr_1: N in structural material in the layer that is
                decomposing
            struce_lyr_2: P in structural material in the layer that is
                decomposing
            rnew_lyr_1_1: required C/N ratio of material decomposing to SOM1
            rnew_lyr_2_1: required C/P ratio of material decomposing to SOM1
            rnew_lyr_1_2: required C/N ratio of material decomposing to SOM2
            rnew_lyr_2_2: required C/P ratio of material decomposing to SOM2
            minerl_1_1: surface mineral N
            minerl_1_2: surface mineral P

        Returns:
            The change in one state variable, where the state variable is
                specified by return_type:
                d_strucc_lyr, change in C in the decomposing layer, if
                    return_type is 'd_strucc'
                d_struce_lyr_1, change in N in the decomposing lyr, if
                    return_type is 'd_struce_1'
                d_struce_lyr_2, change in P in the decomposing lyr, if
                    return_type is 'd_struce_2'
                d_minerl_1_1, change in surface mineral N, if return_type is
                    'd_minerl_1_1'
                d_minerl_1_2, change in surface mineral P, if return_type is
                    'd_minerl_1_2'
                d_gromin_1, change in gross N mineralization, if return_type is
                    'd_gromin_1'
                d_som2c_lyr, change in C in SOM2, if return_type is 'd_som2c'
                d_som2e_lyr_1, change in N in SOM2, if return_type is
                    'd_som2e_1'
                d_som2e_lyr_2, change in P in SOM2, if return_type is
                    'd_som2e_2'
                d_som1c_lyr, change in C in SOM1, if return_type is 'd_som1c'
                d_som1e_lyr_1, change in N in SOM1, if return_type is
                    'd_som1e_1'
                d_som1e_lyr_2, change in P in SOM1, if return_type is
                    'd_som1e_2'
        """
        # initialize change (delta, d) in state variables
        d_strucc_lyr = 0  # change in structural C in decomposing lyr
        d_struce_lyr_1 = 0  # change in structural N in decomposing lyr
        d_struce_lyr_2 = 0  # change in structural P in decomposing lyr
        d_minerl_1_1 = 0  # change in surface mineral N
        d_minerl_1_2 = 0  # change in surface mineral P
        d_gromin_1 = 0  # change in gross N mineralization
        d_som2c_lyr = 0  # change in C in SOM2 in lyr
        d_som2e_lyr_1 = 0  # change in N in SOM2 in lyr
        d_som2e_lyr_2 = 0  # change in P in SOM2 in lyr
        d_som1c_lyr = 0  # change in C in SOM1 in lyr
        d_som1e_lyr_1 = 0  # change in N in SOM1 in lyr
        d_som1e_lyr_2 = 0  # change in P in SOM1 in lyr

        decompose_mask = (
            ((aminrl_1 > 0.0000001) | (
                (strucc_lyr / struce_lyr_1) <= rnew_lyr_1_1)) &
            ((aminrl_2 > 0.0000001) | (
                (strucc_lyr / struce_lyr_2) <= rnew_lyr_2_1)))

        if decompose_mask:
            d_strucc_lyr = -tcflow
            # material decomposes first to som2
            tosom2 = tcflow * ligcon  # line 127 Declig.f

            # respiration associated with decomposition to som2
            co2los = tosom2 * rsplig  # line 130 Declig.f
            mnrflo_1 = co2los * struce_lyr_1 / strucc_lyr  # line 132
            d_struce_lyr_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            if mnrflo_1 > 0:
                d_gromin_1 += mnrflo_1
            mnrflo_2 = co2los * struce_lyr_2 / strucc_lyr
            d_struce_lyr_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            net_tosom2 = tosom2 - co2los  # line 136 Declig.f
            d_som2c_lyr += net_tosom2  # line 140 Declig.f

            # N and P flows from struce_lyr to som2e_lyr, line 145 Declig.f
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom2, strucc_lyr, rnew_lyr_1_2, struce_lyr_1,
                    minerl_1_1)
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom2, strucc_lyr, rnew_lyr_1_2, struce_lyr_1,
                    minerl_1_1)
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom2, strucc_lyr, rnew_lyr_1_2, struce_lyr_1,
                    minerl_1_1)
            # schedule flows
            d_struce_lyr_1 -= material_leaving_a
            d_som2e_lyr_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                d_gromin_1 += mineral_flow

            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom2, strucc_lyr, rnew_lyr_2_2, struce_lyr_2,
                    minerl_1_2)
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom2, strucc_lyr, rnew_lyr_2_2, struce_lyr_2,
                    minerl_1_2)
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom2, strucc_lyr, rnew_lyr_2_2, struce_lyr_2,
                    minerl_1_2)
            # schedule flows
            d_struce_lyr_2 -= material_leaving_a
            d_som2e_lyr_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow

            # what's left decomposes to som1
            tosom1 = tcflow - tosom2  # line 160 Declig.f
            co2los = tosom1 * ps1co2_lyr  # line 163 Declig.f

            # respiration associated with decomposition to som1
            mnrflo_1 = co2los * struce_lyr_1 / strucc_lyr  # line 165
            d_struce_lyr_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            if mnrflo_1 > 0:
                d_gromin_1 += mnrflo_1
            mnrflo_2 = co2los * struce_lyr_2 / strucc_lyr
            d_struce_lyr_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            net_tosom1 = tosom1 - co2los  # line 169 Declig.f
            d_som1c_lyr += net_tosom1  # line 173 Declig.f

            # N and P flows from struce_lyr to som1e_lyr, line 178 Declig.f
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom1, strucc_lyr, rnew_lyr_1_1, struce_lyr_1,
                    minerl_1_1)
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom1, strucc_lyr, rnew_lyr_1_1, struce_lyr_1,
                    minerl_1_1)
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom1, strucc_lyr, rnew_lyr_1_1, struce_lyr_1,
                    minerl_1_1)
            # schedule flows
            d_struce_lyr_1 -= material_leaving_a
            d_som1e_lyr_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                d_gromin_1 += mineral_flow

            # P
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom1, strucc_lyr, rnew_lyr_2_1, struce_lyr_2,
                    minerl_1_2)
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom1, strucc_lyr, rnew_lyr_2_1, struce_lyr_2,
                    minerl_1_2)
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom1, strucc_lyr, rnew_lyr_2_1, struce_lyr_2,
                    minerl_1_2)
            # schedule flows
            d_struce_lyr_2 -= material_leaving_a
            d_som1e_lyr_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow

        if return_type == 'd_strucc':
            return d_strucc_lyr
        elif return_type == 'd_struce_1':
            return d_struce_lyr_1
        elif return_type == 'd_struce_2':
            return d_struce_lyr_2
        elif return_type == 'd_minerl_1_1':
            return d_minerl_1_1
        elif return_type == 'd_minerl_1_2':
            return d_minerl_1_2
        elif return_type == 'd_gromin_1':
            return d_gromin_1
        elif return_type == 'd_som2c':
            return d_som2c_lyr
        elif return_type == 'd_som2e_1':
            return d_som2e_lyr_1
        elif return_type == 'd_som2e_2':
            return d_som2e_lyr_2
        elif return_type == 'd_som1c':
            return d_som1c_lyr
        elif return_type == 'd_som1e_1':
            return d_som1e_lyr_1
        elif return_type == 'd_som1e_2':
            return d_som1e_lyr_2
    return _declig


def decomposition_point(
        inputs, params, state_var, year_reg, month_reg, pp_reg, rnew_dict,
        pevap):
    """Point implementation of decomposition.

    Parameters:
        inputs (dict): dictionary of input values, including precipitation,
            temperature, and soil pH
        params (dict): dictionary of parameter values
        state_var (dict): dictionary of state variables prior to
            decomposition
        year_reg (dict): dictionary of values that are updated once per year,
            including annual precipitation and annual N deposition
        month_reg (dict): dictionary of values that are shared between
            submodels, including snowmelt
        pp_reg (dict): dictionary of persistent parameters, calculated upon
            model initialization
        rnew_dict (dict): dictionary of required ratios for aboveground
            decomposition, calculated during initialization
        evap (float): reference evapotranspiration

    Returns
        dictionary of modified state variables following decomposition
    """
    gromin_1 = 0  # gross mineralization of N, used outside decomposition
    rprpet = rprpet_point(
        pevap, month_reg['snowmelt'], state_var['avh2o_3'], inputs['precip'])
    defac = defac_point(
        state_var['snow'], inputs['min_temp'], inputs['max_temp'], rprpet,
        params['teff_1'], params['teff_2'], params['teff_3'], params['teff_4'])
    anerb = calc_anerb_point(
        rprpet, pevap, params['drain'], params['aneref_1'], params['aneref_2'],
        params['aneref_3'])

    # pH effect on decomposition for structural material, line 45
    pheff_struc = numpy.clip(
        (0.5 + (1.1 / numpy.pi) *
            numpy.arctan(numpy.pi * 0.7 * (inputs['pH'] - 4.))), 0, 1)

    # pH effect on decomposition for metabolic material, line 158
    pheff_metab = numpy.clip(
        (0.5 + (1.14 / numpy.pi) *
            numpy.arctan(numpy.pi * 0.7 * (inputs['pH'] - 4.8))), 0, 1)

    # calculate aminrl_1, intermediate surface mineral N that is tracked
    # during decomposition
    aminrl_1 = state_var['minerl_1_1']

    # calculate aminrl_2, intermediate surface mineral P that is tracked
    # during decomposition
    fsol = fsfunc_point(
        state_var['minerl_1_2'], params['pslsrb'], params['sorpmx'])
    aminrl_2 = state_var['minerl_1_2'] * fsol

    # monthly N fixation
    state_var['minerl_1_1'] = monthly_N_fixation_point(
        inputs['precip'], year_reg['annual_precip'], year_reg['baseNdep'],
        params['epnfs_2'], state_var['minerl_1_1'])

    for _ in xrange(1):  # TODO eventually should be xrange(4)
        # initialize change (delta, d) in state variables for this decomp step
        d_minerl_1_1 = 0  # change in surface mineral N
        d_minerl_1_2 = 0  # change in surface mineral P

        d_strucc_1 = 0  # change in surface structural C
        d_struce_1_1 = 0  # change in surface structural N
        d_struce_1_2 = 0  # change in surface strctural P
        d_som2c_1 = 0  # change in surface SOM2 C
        d_som2e_1_1 = 0  # change in surface SOM2 N
        d_som2e_1_2 = 0  # change in surface SOM2 P
        d_som1c_1 = 0  # change in surface SOM1 C
        d_som1e_1_1 = 0  # change in surface SOM1 N
        d_som1e_1_2 = 0  # change in surface SOM1 P

        d_strucc_2 = 0  # change in soil structural C
        d_struce_2_1 = 0  # change in soil structural N
        d_struce_2_2 = 0  # change in soil structural P
        d_som2c_2 = 0  # change in soil SOM2 C
        d_som2e_2_1 = 0  # change in soil SOM2 N
        d_som2e_2_2 = 0  # change in soil SOM2 P
        d_som1c_2 = 0  # change in soil SOM1 C
        d_som1e_2_1 = 0  # change in soil SOM1 N
        d_som1e_2_2 = 0  # change in soil SOM1 P

        d_metabc_1 = 0  # change in surface metabolic C
        d_metabe_1_1 = 0  # change in surface metabolic N
        d_metabe_1_2 = 0  # change in surface metabolic P

        d_metabc_2 = 0  # change in soil metabolic C
        d_metabe_2_1 = 0  # change in soil metabolic N
        d_metabe_2_2 = 0  # change in soil metabolic P

        d_som3c = 0  # change in passive organic C
        d_som3e_1 = 0  # change in passive organic N
        d_som3e_2 = 0  # change in passive organic P

        d_parent_2 = 0  # change in parent P
        d_secndy_2 = 0  # change in secondary P
        d_occlud = 0  # change in occluded P

        d_minerl_P_dict = {}
        for lyr in xrange(1, params['nlayer'] + 1):
            d_minerl_P_dict['d_minerl_{}_2'.format(lyr)] = 0

        # litdec.f
        # decomposition of strucc_1, surface structural material
        tcflow = (min(
            state_var['strucc_1'], params['strmax_1']) * defac *
            params['dec1_1'] *
            math.exp(-params['pligst_1'] * state_var['strlig_1']) *
            0.020833 * pheff_struc)

        d_strucc_1 += declig_point(
            'd_strucc')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_struce_1_1 += declig_point(
            'd_struce_1')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_struce_1_2 += declig_point(
            'd_struce_2')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_minerl_1_1 += declig_point(
            'd_minerl_1_1')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_minerl_1_2 += declig_point(
            'd_minerl_1_2')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        gromin_1 += declig_point(
            'd_gromin_1')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_som2c_1 += declig_point(
            'd_som2c')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_som2e_1_1 += declig_point(
            'd_som2e_1')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_som2e_1_2 += declig_point(
            'd_som2e_2')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_som1c_1 += declig_point(
            'd_som1c')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_som1e_1_1 += declig_point(
            'd_som1e_1')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])
        d_som1e_1_2 += declig_point(
            'd_som1e_2')(
                aminrl_1, aminrl_2, state_var['strlig_1'], params['rsplig'],
                params['ps1co2_1'], state_var['strucc_1'],
                tcflow, state_var['struce_1_1'],
                state_var['struce_1_2'], rnew_dict['rnewas_1_1'],
                rnew_dict['rnewas_2_1'], rnew_dict['rnewas_1_2'],
                rnew_dict['rnewas_2_2'], state_var['minerl_1_1'],
                state_var['minerl_1_2'])

        # decomposition of strucc_2, soil structural material: line 99 Litdec.f
        tcflow = (
            min(state_var['strucc_2'], params['strmax_2']) * defac *
            params['dec1_2'] *
            math.exp(-params['pligst_2'] * state_var['strlig_2']) *
            anerb * 0.020833 * pheff_struc)

        d_strucc_2 += declig_point(
            'd_strucc')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_struce_2_1 += declig_point(
            'd_struce_1')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_struce_2_2 += declig_point(
            'd_struce_2')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_minerl_1_1 += declig_point(
            'd_minerl_1_1')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_minerl_1_2 += declig_point(
            'd_minerl_1_2')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        gromin_1 += declig_point(
            'd_gromin_1')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_som2c_2 += declig_point(
            'd_som2c')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_som2e_2_1 += declig_point(
            'd_som2e_1')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_som2e_2_2 += declig_point(
            'd_som2e_2')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_som1c_2 += declig_point(
            'd_som1c')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_som1e_2_1 += declig_point(
            'd_som1e_1')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])
        d_som1e_2_2 += declig_point(
            'd_som1e_2')(
                aminrl_1, aminrl_2, state_var['strlig_2'], params['rsplig'],
                params['ps1co2_2'], state_var['strucc_2'], tcflow,
                state_var['struce_2_1'], state_var['struce_2_2'],
                rnew_dict['rnewbs_1_1'], rnew_dict['rnewbs_2_1'],
                rnew_dict['rnewbs_1_2'], rnew_dict['rnewbs_2_2'],
                state_var['minerl_1_1'], state_var['minerl_1_2'])

        # decomposition of surface metabolic material: line 136 Litdec.f
        # C/N ratio for surface metabolic residue
        rceto1_1 = agdrat_point(
            state_var['metabe_1_1'], state_var['metabc_1'],
            params['pcemic1_1_1'], params['pcemic1_2_1'],
            params['pcemic1_3_1'])
        # C/P ratio for surface metabolic residue
        rceto1_2 = agdrat_point(
            state_var['metabe_1_2'], state_var['metabc_1'],
            params['pcemic1_1_2'], params['pcemic1_2_2'],
            params['pcemic1_3_2'])
        decompose_mask = (
            ((aminrl_1 > 0.0000001) | (
                (state_var['metabc_1'] / state_var['metabe_1_1']) <=
                rceto1_1)) &
            ((aminrl_2 > 0.0000001) | (
                (state_var['metabc_1'] / state_var['metabe_1_2']) <=
                rceto1_2)))  # line 194 Litdec.f
        if decompose_mask:
            tcflow = numpy.clip(
                (state_var['metabc_1'] * defac * params['dec2_1'] * 0.020833 *
                    pheff_metab), 0,
                state_var['metabc_1'])
            co2los = tcflow * params['pmco2_1']
            d_metabc_1 -= tcflow
            # respiration, line 201 Litdec.f
            mnrflo_1 = (
                co2los * state_var['metabe_1_1'] / state_var['metabc_1'])
            d_metabe_1_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            gromin_1 += mnrflo_1
            mnrflo_2 = (
                co2los * state_var['metabe_1_2'] / state_var['metabc_1'])
            d_metabe_1_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            net_tosom1 = tcflow - co2los  # line 210 Litdec.f
            d_som1c_1 += net_tosom1
            # N and P flows from metabe_1 to som1e_1, line 222 Litdec.f
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom1, state_var['metabc_1'], rceto1_1,
                    state_var['metabe_1_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom1, state_var['metabc_1'], rceto1_1,
                    state_var['metabe_1_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom1, state_var['metabc_1'], rceto1_1,
                    state_var['metabe_1_1'], state_var['minerl_1_1'])
            # schedule flows
            d_metabe_1_1 -= material_leaving_a
            d_som1e_1_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow

            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom1, state_var['metabc_1'], rceto1_2,
                    state_var['metabe_1_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom1, state_var['metabc_1'], rceto1_2,
                    state_var['metabe_1_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom1, state_var['metabc_1'], rceto1_2,
                    state_var['metabe_1_2'], state_var['minerl_1_2'])
            # schedule flows
            d_metabe_1_2 -= material_leaving_a
            d_som1e_1_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow

        # decomposition of soil metabolic material: line 136 Litdec.f
        # C/N ratio for belowground material to som1
        rceto1_1 = bgdrat_point(
            aminrl_1, params['varat1_1_1'], params['varat1_2_1'],
            params['varat1_3_1'])
        # C/P ratio for soil metabolic material
        rceto1_2 = bgdrat_point(
            aminrl_2, params['varat1_1_2'], params['varat1_2_2'],
            params['varat1_3_2'])
        decompose_mask = (
            ((aminrl_1 > 0.0000001) | (
                (state_var['metabc_2'] / state_var['metabe_2_1']) <=
                rceto1_1)) &
            ((aminrl_2 > 0.0000001) | (
                (state_var['metabc_2'] / state_var['metabe_2_2']) <=
                rceto1_2)))  # line 194 Litdec.f
        if decompose_mask:
            tcflow = numpy.clip(
                (state_var['metabc_2'] * defac * params['dec2_2'] * 0.020833 *
                    pheff_metab * anerb),
                0, state_var['metabc_2'])
            co2los = tcflow * params['pmco2_2']
            d_metabc_2 -= tcflow
            # respiration, line 201 Litdec.f
            mnrflo_1 = co2los * state_var['metabe_2_1'] / state_var['metabc_2']
            d_metabe_2_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            gromin_1 += mnrflo_1
            mnrflo_2 = co2los * state_var['metabe_2_2'] / state_var['metabc_2']
            d_metabe_2_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            net_tosom1 = tcflow - co2los  # line 210 Litdec.f
            d_som1c_2 += net_tosom1
            # N and P flows from metabe_2 to som1e_2, line 222 Litdec.f
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom1, state_var['metabc_2'], rceto1_1,
                    state_var['metabe_2_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom1, state_var['metabc_2'], rceto1_1,
                    state_var['metabe_2_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom1, state_var['metabc_2'], rceto1_1,
                    state_var['metabe_2_1'], state_var['minerl_1_1'])
            # schedule flows
            d_metabe_2_1 -= material_leaving_a
            d_som1e_2_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow

            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom1, state_var['metabc_2'], rceto1_2,
                    state_var['metabe_2_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom1, state_var['metabc_2'], rceto1_2,
                    state_var['metabe_2_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom1, state_var['metabc_2'], rceto1_2,
                    state_var['metabe_2_2'], state_var['minerl_1_2'])
            # schedule flows
            d_metabe_2_2 -= material_leaving_a
            d_som1e_2_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow

        # somdec.f
        # surface SOM1 to surface SOM2
        # rceto2_iel: required ratio for flow from surface SOM1 to SOM2
        radds1_1 = (
            params['rad1p_1_1'] + params['rad1p_2_1'] *
            ((state_var['som1c_1'] / state_var['som1e_1_1']) -
                params['pcemic1_2_1']))
        rceto2_1_surface = max(
            (state_var['som1c_1'] / state_var['som1e_1_1'] + radds1_1),
            params['rad1p_3_1'])

        radds1_2 = (
            params['rad1p_1_2'] + params['rad1p_2_2'] *
            ((state_var['som1c_1'] / state_var['som1e_1_2']) -
                params['pcemic1_2_2']))
        rceto2_2_surface = max(
            (state_var['som1c_1'] / state_var['som1e_1_2'] + radds1_2),
            params['rad1p_3_2'])

        decompose_mask = (
            ((aminrl_1 > 0.0000001) | (
                (state_var['som1c_1'] / state_var['som1e_1_1']) <=
                rceto2_1_surface)) &
            ((aminrl_2 > 0.0000001) | (
                (state_var['som1c_1'] / state_var['som1e_1_2']) <=
                rceto2_2_surface)))  # line 92
        if decompose_mask:
            tcflow = (
                state_var['som1c_1'] * defac * params['dec3_1'] * 0.020833 *
                pheff_struc)
            co2los = tcflow * params['p1co2a_1']
            d_som1c_1 -= tcflow
            # respiration, line 105 Somdec.f
            mnrflo_1 = co2los * state_var['som1e_1_1'] / state_var['som1c_1']
            d_som1e_1_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            gromin_1 += mnrflo_1
            mnrflo_2 = co2los * state_var['som1e_1_2'] / state_var['som1c_1']
            d_som1e_1_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            net_tosom2 = tcflow - co2los
            d_som2c_1 += net_tosom2
            # N and P flows from som1e_1 to som2e_1, line 123 Somdec.f
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom2, state_var['som1c_1'], rceto2_1_surface,
                    state_var['som1e_1_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom2, state_var['som1c_1'], rceto2_1_surface,
                    state_var['som1e_1_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom2, state_var['som1c_1'], rceto2_1_surface,
                    state_var['som1e_1_1'], state_var['minerl_1_1'])
            # schedule flows
            d_som1e_1_1 -= material_leaving_a
            d_som2e_1_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow

            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom2, state_var['som1c_1'], rceto2_2_surface,
                    state_var['som1e_1_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom2, state_var['som1c_1'], rceto2_2_surface,
                    state_var['som1e_1_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom2, state_var['som1c_1'], rceto2_2_surface,
                    state_var['som1e_1_2'], state_var['minerl_1_2'])
            # schedule flows
            d_som1e_1_2 -= material_leaving_a
            d_som2e_1_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow

        # soil SOM1 to soil SOM3 and SOM2
        # required ratios are those pertaining to decomposition to SOM2
        rceto2_1 = bgdrat_point(
            aminrl_1, params['varat22_1_1'], params['varat22_2_1'],
            params['varat22_3_1'])  # line 141 Somdec.f
        rceto2_2 = bgdrat_point(
            aminrl_2, params['varat22_1_2'], params['varat22_2_2'],
            params['varat22_3_2'])

        decompose_mask = (
            ((aminrl_1 > 0.0000001) | (
                (state_var['som1c_2'] / state_var['som1e_2_1']) <=
                rceto2_1)) &
            ((aminrl_2 > 0.0000001) | (
                (state_var['som1c_2'] / state_var['som1e_2_2']) <=
                rceto2_2)))  # line 171
        if decompose_mask:
            tcflow = (
                state_var['som1c_2'] * defac * params['dec3_2'] *
                pp_reg['eftext'] * anerb * 0.020833 * pheff_metab)
            co2los = tcflow * params['p1co2_2']
            d_som1c_2 -= tcflow
            # respiration, line 179 Somdec.f
            mnrflo_1 = co2los * state_var['som1e_2_1'] / state_var['som1c_2']
            d_som1e_2_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            gromin_1 += mnrflo_1
            mnrflo_2 = co2los * state_var['som1e_2_2'] / state_var['som1c_2']
            d_som1e_2_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            tosom3 = (
                tcflow * pp_reg['fps1s3'] *
                (1. + params['animpt'] * (1. - anerb)))
            d_som3c += tosom3
            # C/<iel> ratios of material entering som3e
            rceto3_1 = bgdrat_point(
                aminrl_1, params['varat3_1_1'], params['varat3_2_1'],
                params['varat3_3_1'])
            rceto3_2 = bgdrat_point(
                aminrl_2, params['varat3_1_2'], params['varat3_2_2'],
                params['varat3_3_2'])
            # N and P flows from soil som1e to som3e, line 198
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom3, state_var['som1c_2'], rceto3_1,
                    state_var['som1e_2_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom3, state_var['som1c_2'], rceto3_1,
                    state_var['som1e_2_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom3, state_var['som1c_2'], rceto3_1,
                    state_var['som1e_2_1'], state_var['minerl_1_1'])
            # schedule flows
            d_som1e_2_1 -= material_leaving_a
            d_som3e_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow
            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom3, state_var['som1c_2'], rceto3_2,
                    state_var['som1e_2_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom3, state_var['som1c_2'], rceto3_2,
                    state_var['som1e_2_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom3, state_var['som1c_2'], rceto3_2,
                    state_var['som1e_2_2'], state_var['minerl_1_2'])
            # schedule flows
            d_som1e_2_2 -= material_leaving_a
            d_som3e_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow

            # organic leaching: line 204 Somdec.f
            if month_reg['amov_2'] > 0:
                linten = min(
                    (1. - (params['omlech_3'] - month_reg['amov_2']) /
                        params['omlech_3']), 1.)
                cleach = tcflow * pp_reg['orglch'] * linten
                # N leaching: line 230
                rceof1_1 = (state_var['som1c_2'] / state_var['som1e_2_1']) * 2
                orgflow_1 = cleach / rceof1_1
                d_som1e_2_1 -= orgflow_1
                # P leaching: line 232
                rceof1_2 = (state_var['som1c_2'] / state_var['som1e_2_2']) * 35
                orgflow_2 = cleach / rceof1_2
                d_som1e_2_2 -= orgflow_2
            else:
                cleach = 0

            net_tosom2 = tcflow - co2los - tosom3 - cleach
            d_som2c_2 += net_tosom2
            # N and P flows from som1e_2 to som2e_2, line 257
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom2, state_var['som1c_2'], rceto2_1,
                    state_var['som1e_2_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom2, state_var['som1c_2'], rceto2_1,
                    state_var['som1e_2_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom2, state_var['som1c_2'], rceto2_1,
                    state_var['som1e_2_1'], state_var['minerl_1_1'])
            # schedule flows
            d_som1e_2_1 -= material_leaving_a
            d_som2e_2_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow
            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    net_tosom2, state_var['som1c_2'], rceto2_2,
                    state_var['som1e_2_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    net_tosom2, state_var['som1c_2'], rceto2_2,
                    state_var['som1e_2_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    net_tosom2, state_var['som1c_2'], rceto2_2,
                    state_var['som1e_2_2'], state_var['minerl_1_2'])
            # schedule flows
            d_som1e_2_2 -= material_leaving_a
            d_som2e_2_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow
        # Soil SOM2 decomposing to soil SOM1 and SOM3, line 269 Somdec.f
        decompose_mask = (
            ((aminrl_1 > 0.0000001) | (
                (state_var['som2c_2'] / state_var['som2e_2_1']) <=
                rceto1_1)) &
            ((aminrl_2 > 0.0000001) | (
                (state_var['som2c_2'] / state_var['som2e_2_2']) <=
                rceto1_2)))  # line 298
        if decompose_mask:
            tcflow = (
                state_var['som2c_2'] * defac * params['dec5_2'] * anerb *
                0.020833 * pheff_metab)
            co2los = tcflow * params['p2co2_2']
            d_som2c_2 -= tcflow
            # respiration, line 304 Somdec.f
            mnrflo_1 = co2los * state_var['som2e_2_1'] / state_var['som2c_2']
            d_som2e_2_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            gromin_1 += mnrflo_1
            mnrflo_2 = co2los * state_var['som2e_2_2'] / state_var['som2c_2']
            d_som2e_2_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            tosom3 = (
                tcflow * pp_reg['fps2s3'] *
                (1. + params['animpt'] * (1.0 - anerb)))
            d_som3c += tosom3
            # N and P flows from soil som2e to som3e
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom3, state_var['som2c_2'], rceto3_1,
                    state_var['som2e_2_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom3, state_var['som2c_2'], rceto3_1,
                    state_var['som2e_2_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom3, state_var['som2c_2'], rceto3_1,
                    state_var['som2e_2_1'], state_var['minerl_1_1'])
            # schedule flows
            d_som2e_2_1 -= material_leaving_a
            d_som3e_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow
            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom3, state_var['som2c_2'], rceto3_2,
                    state_var['som2e_2_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom3, state_var['som2c_2'], rceto3_2,
                    state_var['som2e_2_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom3, state_var['som2c_2'], rceto3_2,
                    state_var['som2e_2_2'], state_var['minerl_1_2'])
            # schedule flows
            d_som2e_2_2 -= material_leaving_a
            d_som3e_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow

            # rest of the flow from SOM2 goes to SOM1
            tosom1 = tcflow - co2los - tosom3  # line 333 Somdec.f
            d_som1c_2 += tosom1
            # N and P flows from som2e_2 to som1e_2, line 344
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom1, state_var['som2c_2'], rceto1_1,
                    state_var['som2e_2_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom1, state_var['som2c_2'], rceto1_1,
                    state_var['som2e_2_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom1, state_var['som2c_2'], rceto1_1,
                    state_var['som2e_2_1'], state_var['minerl_1_1'])
            # schedule flows
            d_som2e_2_1 -= material_leaving_a
            d_som1e_2_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow
            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom1, state_var['som2c_2'], rceto1_2,
                    state_var['som2e_2_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom1, state_var['som2c_2'], rceto1_2,
                    state_var['som2e_2_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom1, state_var['som2c_2'], rceto1_2,
                    state_var['som2e_2_2'], state_var['minerl_1_2'])
            # schedule flows
            d_som2e_2_2 -= material_leaving_a
            d_som1e_2_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow
        # Surface SOM2 decomposes to surface SOM1
        # ratios of material decomposing to SOM1
        decompose_mask = (
            ((aminrl_1 > 0.0000001) | (
                (state_var['som2c_1'] / state_var['som2e_1_1']) <=
                rceto1_1)) &
            ((aminrl_2 > 0.0000001) | (
                (state_var['som2c_1'] / state_var['som2e_1_2']) <=
                rceto1_2)))  # line 171
        if decompose_mask:
            tcflow = (
                state_var['som2c_1'] * defac * params['dec5_1'] * 0.020833 *
                pheff_struc)
            co2los = tcflow * params['p2co2_1']  # line 385
            d_som2c_1 -= tcflow
            # respiration, line 388
            mnrflo_1 = co2los * state_var['som2e_1_1'] / state_var['som2c_1']
            d_som2e_1_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            gromin_1 += mnrflo_1
            mnrflo_2 = co2los * state_var['som2e_1_2'] / state_var['som2c_1']
            d_som2e_1_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            tosom1 = tcflow - co2los  # line 393
            d_som1c_1 += tosom1
            # N and P flows from surface som2e to surface som1e, line 404
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom1, state_var['som2c_1'], rceto1_1,
                    state_var['som2e_1_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom1, state_var['som2c_1'], rceto1_1,
                    state_var['som2e_1_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom1, state_var['som2c_1'], rceto1_1,
                    state_var['som2e_1_1'], state_var['minerl_1_1'])
            # schedule flows
            d_som2e_1_1 -= material_leaving_a
            d_som1e_1_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow
            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom1, state_var['som2c_1'], rceto1_2,
                    state_var['som2e_1_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom1, state_var['som2c_1'], rceto1_2,
                    state_var['som2e_1_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom1, state_var['som2c_1'], rceto1_2,
                    state_var['som2e_1_2'], state_var['minerl_1_2'])
            # schedule flows
            d_som2e_1_2 -= material_leaving_a
            d_som1e_1_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow
        # SOM3 decomposes to soil SOM1
        pheff_som3 = numpy.clip(
            (0.5 + (1.1 / numpy.pi) *
                numpy.arctan(numpy.pi * 0.7 * (inputs['pH'] - 3.))), 0, 1)
        decompose_mask = (
            ((aminrl_1 > 0.0000001) | (
                (state_var['som3c'] / state_var['som3e_1']) <= rceto1_1)) &
            ((aminrl_2 > 0.0000001) | (
                (state_var['som3c'] / state_var['som3e_2']) <= rceto1_2)))
        if decompose_mask:
            tcflow = (
                state_var['som3c'] * defac * params['dec4'] * anerb *
                0.020833 * pheff_som3)
            co2los = tcflow * params['p3co2'] * anerb  # line 442
            d_som3c -= tcflow
            # respiration, line 446
            mnrflo_1 = co2los * state_var['som3e_1'] / state_var['som3c']
            d_som3e_1 -= mnrflo_1
            d_minerl_1_1 += mnrflo_1
            gromin_1 += mnrflo_1
            mnrflo_2 = co2los * state_var['som3e_2'] / state_var['som3c']
            d_som3e_2 -= mnrflo_2
            d_minerl_1_2 += mnrflo_2

            tosom1 = tcflow - co2los
            d_som1c_2 += tosom1
            # N and P flows from som3e to som1e_2, line 461 Somdec.f
            # N first
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom1, state_var['som3c'], rceto1_1,
                    state_var['som3e_1'], state_var['minerl_1_1'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom1, state_var['som3c'], rceto1_1,
                    state_var['som3e_1'], state_var['minerl_1_1'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom1, state_var['som3c'], rceto1_1,
                    state_var['som3e_1'], state_var['minerl_1_1'])
            # schedule flows
            d_som3e_1 -= material_leaving_a
            d_som1e_2_1 += material_arriving_b
            d_minerl_1_1 += mineral_flow
            if mineral_flow > 0:
                gromin_1 += mineral_flow
            # P second
            material_leaving_a = esched_point(
                'material_leaving_a')(
                    tosom1, state_var['som3c'], rceto1_2,
                    state_var['som3e_2'], state_var['minerl_1_2'])
            material_arriving_b = esched_point(
                'material_arriving_b')(
                    tosom1, state_var['som3c'], rceto1_2,
                    state_var['som3e_2'], state_var['minerl_1_2'])
            mineral_flow = esched_point(
                'mineral_flow')(
                    tosom1, state_var['som3c'], rceto1_2,
                    state_var['som3e_2'], state_var['minerl_1_2'])
            # schedule flows
            d_som3e_2 -= material_leaving_a
            d_som1e_2_2 += material_arriving_b
            d_minerl_1_2 += mineral_flow
        # Surface SOM2 flows to soil SOM2 via mixing
        tcflow = state_var['som2c_1'] * params['cmix'] * defac * 0.020833
        d_som2c_1 -= tcflow
        d_som2c_2 += tcflow
        # N and P flows from som2e_1 to som2e_2, line 495 Somdec.f
        # N first
        mix_ratio_1 = state_var['som2c_1'] / state_var['som2e_1_1']  # line 495
        material_leaving_a = esched_point(
            'material_leaving_a')(
                tcflow, state_var['som2c_1'], mix_ratio_1,
                state_var['som2e_1_1'], state_var['minerl_1_1'])
        material_arriving_b = esched_point(
            'material_arriving_b')(
                tcflow, state_var['som2c_1'], mix_ratio_1,
                state_var['som2e_1_1'], state_var['minerl_1_1'])
        mineral_flow = esched_point(
            'mineral_flow')(
                tcflow, state_var['som2c_1'], mix_ratio_1,
                state_var['som2e_1_1'], state_var['minerl_1_1'])
        # schedule flows
        d_som2e_1_1 -= material_leaving_a
        d_som2e_2_1 += material_arriving_b
        d_minerl_1_1 += mineral_flow
        # P second
        mix_ratio_2 = state_var['som2c_1'] / state_var['som2e_1_2']  # line 495
        material_leaving_a = esched_point(
            'material_leaving_a')(
                tcflow, state_var['som2c_1'], mix_ratio_2,
                state_var['som2e_1_2'], state_var['minerl_1_2'])
        material_arriving_b = esched_point(
            'material_arriving_b')(
                tcflow, state_var['som2c_1'], mix_ratio_2,
                state_var['som2e_1_2'], state_var['minerl_1_2'])
        mineral_flow = esched_point(
            'mineral_flow')(
                tcflow, state_var['som2c_1'], mix_ratio_2,
                state_var['som2e_1_2'], state_var['minerl_1_2'])
        # schedule flows
        d_som2e_1_2 -= material_leaving_a
        d_som2e_2_2 += material_arriving_b
        d_minerl_1_2 += mineral_flow

        # P mineral flows: Pschem.f
        # flow from parent to mineral: line 141
        fparnt = params['pparmn_2'] * state_var['parent_2'] * defac * 0.020833
        d_minerl_1_2 += fparnt
        d_parent_2 -= fparnt

        # flow from secondary to mineral: line 158
        fsecnd = params['psecmn_2'] * state_var['secndy_2'] * defac * 0.020833
        d_minerl_1_2 += fsecnd
        d_secndy_2 -= fsecnd

        # flow from mineral to secondary: line 163
        fmnsec = (
            params['pmnsec_2'] * state_var['minerl_1_2'] * (1 - fsol) * defac *
            0.020833)
        d_minerl_1_2 -= fmnsec
        d_secndy_2 += fmnsec
        for lyr in xrange(2, params['nlayer'] + 1):
            fmnsec = (
                params['pmnsec_2'] * state_var['minerl_{}_2'.format(lyr)] *
                (1 - fsol) * defac * 0.020833)
            d_minerl_P_dict['d_minerl_{}_2'.format(lyr)] -= fmnsec
            d_secndy_2 += fmnsec

        # flow from secondary to occluded: line 171
        fsecoc = params['psecoc1'] * state_var['secndy_2'] * defac * 0.020833
        d_secndy_2 -= fsecoc
        d_occlud += fsecoc

        # flow from occluded to secondary
        focsec = params['psecoc2'] * state_var['occlud'] * defac * 0.020833
        d_occlud -= focsec
        d_secndy_2 += focsec

        # update state variables: perform flows calculated in previous lines
        state_var['minerl_1_1'] += d_minerl_1_1
        state_var['minerl_1_2'] += d_minerl_1_2

        state_var['strucc_1'] += d_strucc_1
        state_var['struce_1_1'] += d_struce_1_1
        state_var['struce_1_2'] += d_struce_1_2
        state_var['som2c_1'] += d_som2c_1
        state_var['som2e_1_1'] += d_som2e_1_1
        state_var['som2e_1_2'] += d_som2e_1_2
        state_var['som1c_1'] += d_som1c_1
        state_var['som1e_1_1'] += d_som1e_1_1
        state_var['som1e_1_2'] += d_som1e_1_2

        state_var['strucc_2'] += d_strucc_2
        state_var['struce_2_1'] += d_struce_2_1
        state_var['struce_2_2'] += d_struce_2_2
        state_var['som2c_2'] += d_som2c_2
        state_var['som2e_2_1'] += d_som2e_2_1
        state_var['som2e_2_2'] += d_som2e_2_2
        state_var['som1c_2'] += d_som1c_2
        state_var['som1e_2_1'] += d_som1e_2_1
        state_var['som1e_2_2'] += d_som1e_2_2

        state_var['metabc_1'] += d_metabc_1
        state_var['metabe_1_1'] += d_metabe_1_1
        state_var['metabe_1_2'] += d_metabe_1_2

        state_var['metabc_2'] += d_metabc_2
        state_var['metabe_2_1'] += d_metabe_2_1
        state_var['metabe_2_2'] += d_metabe_2_2

        state_var['parent_2'] += d_parent_2
        state_var['secndy_2'] += d_secndy_2
        state_var['occlud'] += d_occlud
        for lyr in xrange(2, params['nlayer'] + 1):
            state_var['minerl_{}_2'.format(lyr)] += (
                d_minerl_P_dict['d_minerl_{}_2'.format(lyr)])

        # update aminrl_1
        aminrl_1 = aminrl_1 + state_var['minerl_1_1'] / 2.

        # update aminrl_2
        fsol = fsfunc_point(
            state_var['minerl_1_2'], params['pslsrb'], params['sorpmx'])
        aminrl_2 = aminrl_2 + (state_var['minerl_1_2'] * fsol) / 2.

    # Calculate volatilization loss of nitrogen as a function of
    # gross mineralization: line 323 Simsom.f
    # volgm = params['vlossg'] * gromin_1
    # state_var['minerl_1_1'] -= volgm
    return state_var


def agdrat_point(anps, tca, pcemic_1_iel, pcemic_2_iel, pcemic_3_iel):
    """Point implementation of `Agdrat.f`.

    Calculate the C/<iel> ratio of new material that is the result of
    decomposition into "box B".

    Parameters:
        anps: <iel> (N or P) in the decomposing stock
        tca: total C in the decomposing stock
        pcemic_1_iel: maximum C/<iel> of new SOM1
        pcemic_2_iel: minimum C/<iel> of new SOM1
        pcemic_3_iel: minimum <iel> content of decomposing material that gives
            minimum C/<iel> of new material

    Returns:
        agdrat, the C/<iel> ratio of new material
    """
    cemicb = (pcemic_2_iel - pcemic_1_iel) / pcemic_3_iel
    if ((tca * 2.5) <= 0.0000000001):
        econt = 0
    else:
        econt = anps / (tca * 2.5)
    if econt > pcemic_3_iel:
        agdrat = pcemic_2_iel
    else:
        agdrat = pcemic_1_iel + econt * cemicb
    return agdrat


def fsfunc_point(minerl_1_2, pslsrb, sorpmx):
        """Calculate the fraction of mineral P that is in solution.

        The fraction of P in solution is influenced by two soil properties:
        the maximum sorption potential of the soil and sorption affinity.

        Parameters:
            minerl_1_2 (float): state variable, surface mineral P
            pslsrb (float): parameter, P sorption affinity
            sorpmx (float): parameter, maximum P sorption of the soil

        Returns:
            fsol, fraction of P in solution
        """
        if minerl_1_2 == 0:
            return 0
        c = sorpmx * (2. - pslsrb) / 2.
        b = sorpmx - minerl_1_2 + c
        labile = (-b + numpy.sqrt(b * b + 4 * c * minerl_1_2)) / 2.
        fsol = labile / minerl_1_2
        return fsol


class foragetests(unittest.TestCase):
    """Regression tests for InVEST forage model."""

    def setUp(self):
        """Create temporary workspace directory."""
        self.workspace_dir = tempfile.mkdtemp()
        self.PROCESSING_DIR = os.path.join(
            self.workspace_dir, "temporary_files")

    def tearDown(self):
        """Clean up remaining files."""
        shutil.rmtree(self.workspace_dir)

    @staticmethod
    def generate_base_args(workspace_dir):
        """Generate a base sample args dict for forage model."""
        args = {
            'workspace_dir': workspace_dir,
            'results_suffix': "",
            'starting_month': 1,
            'starting_year': 2016,
            'n_months': 22,
            'aoi_path': os.path.join(
                SAMPLE_DATA, 'Manlai_soum_WGS84.shp'),
            'bulk_density_path': os.path.join(
                SAMPLE_DATA, 'soil', 'bulkd.tif'),
            'ph_path': os.path.join(
                SAMPLE_DATA, 'soil', 'phihox_sl3.tif'),
            'clay_proportion_path': os.path.join(
                SAMPLE_DATA, 'soil', 'clay.tif'),
            'silt_proportion_path': os.path.join(
                SAMPLE_DATA, 'soil', 'silt.tif'),
            'sand_proportion_path': os.path.join(
                SAMPLE_DATA, 'soil', 'sand.tif'),
            'monthly_precip_path_pattern': os.path.join(
                SAMPLE_DATA, 'CHIRPS_div_by_10',
                'chirps-v2.0.<year>.<month>.tif'),
            'min_temp_path_pattern': os.path.join(
                SAMPLE_DATA, 'temp', 'wc2.0_30s_tmin_<month>.tif'),
            'max_temp_path_pattern': os.path.join(
                SAMPLE_DATA, 'temp', 'wc2.0_30s_tmax_<month>.tif'),
            'site_param_table': os.path.join(
                SAMPLE_DATA, 'site_parameters.csv'),
            'site_param_spatial_index_path': os.path.join(
                SAMPLE_DATA, 'site_index.tif'),
            'veg_trait_path': os.path.join(SAMPLE_DATA, 'pft_trait.csv'),
            'veg_spatial_composition_path_pattern': os.path.join(
                SAMPLE_DATA, 'pft<PFT>.tif'),
            'animal_trait_path': os.path.join(
                SAMPLE_DATA, 'animal_trait_table.csv'),
            'animal_mgmt_layer_path': os.path.join(
                SAMPLE_DATA,
                'sheep_units_density_2016_monitoring_area_WGS84.shp'),
            'initial_conditions_dir': os.path.join(
                SAMPLE_DATA, 'initialization_data'),
        }
        return args

    def assert_all_values_in_raster_within_range(
            self, raster_to_test, minimum_acceptable_value,
            maximum_acceptable_value, nodata_value):
        """Test that `raster_to_test` contains values within acceptable range.

        The values within `raster_to_test` that are not null must be
        greater than or equal to `minimum_acceptable_value` and
        less than or equal to `maximum_acceptable_value`.

        Raises:
            AssertionError if values are outside acceptable range

        Returns:
            None
        """
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                raster_to_test):
            if len(raster_block[raster_block != nodata_value]) == 0:
                continue
            min_val = numpy.amin(
                raster_block[raster_block != nodata_value])
            self.assertGreaterEqual(
                min_val, minimum_acceptable_value,
                msg="Raster contains values smaller than acceptable "
                + "minimum: {}, {} (acceptable min: {})".format(
                    raster_to_test, min_val, minimum_acceptable_value))
            max_val = numpy.amax(
                raster_block[raster_block != nodata_value])
            self.assertLessEqual(
                max_val, maximum_acceptable_value,
                msg="Raster contains values larger than acceptable "
                + "maximum: {}, {} (acceptable max: {})".format(
                    raster_to_test, max_val, maximum_acceptable_value))

    def assert_all_values_in_array_within_range(
            self, array_to_test, minimum_acceptable_value,
            maximum_acceptable_value, nodata_value):
        """Test that `array_to_test` contains values within acceptable range.

        The values within `array_to_test` that are not null must be
        greater than or equal to `minimum_acceptable_value` and
        less than or equal to `maximum_acceptable_value`.

        Raises:
            AssertionError if values are outside acceptable range

        Returns:
            None
        """
        if len(array_to_test[array_to_test != nodata_value]) == 0:
            return
        min_val = numpy.amin(
            array_to_test[array_to_test != nodata_value])
        self.assertGreaterEqual(
            min_val, minimum_acceptable_value,
            msg="Array contains values smaller than acceptable minimum: " +
            "min value: {}, acceptable min: {}".format(
                min_val, minimum_acceptable_value))
        max_val = numpy.amax(
            array_to_test[array_to_test != nodata_value])
        self.assertLessEqual(
            max_val, maximum_acceptable_value,
            msg="Array contains values larger than acceptable maximum: " +
            "max value: {}, acceptable max: {}".format(
                max_val, maximum_acceptable_value))

    @unittest.skip("did not run the whole model, running unit tests only")
    def test_model_runs(self):
        """Test forage model."""
        from natcap.invest import forage

        if not os.path.exists(SAMPLE_DATA):
            self.fail(
                "Sample input directory not found at %s" % SAMPLE_DATA)

        args = foragetests.generate_base_args(self.workspace_dir)
        forage.execute(args)

    def test_shortwave_radiation(self):
        """Test `_shortwave radiation`.

        Use one set of inputs, including a raster with known latitude, to test
        the function `_shortwave radiation` against a result calculated by
        hand.

        Raises:
            AssertionError if the raster created by `_shortwave_radiation`
                contains more than one unique value
            AssertionError if the value returned by `_shortwave_radiation' is
                not within 0.01 of the value calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        fill_value = 0
        template_raster = os.path.join(
            self.workspace_dir, 'template_raster.tif')

        create_constant_raster(template_raster, fill_value)

        month = 5
        shwave_path = os.path.join(self.workspace_dir, 'shwave.tif')
        forage._shortwave_radiation(template_raster, month, shwave_path)

        # assert the value in the raster `shwave_path` is equal to value
        # calculated by hand
        result_set = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                shwave_path):
            result_set.update(numpy.unique(raster_block))
        self.assertEqual(
            len(result_set), 1,
            msg="One unique value expected in shortwave radiation raster")
        test_result = list(result_set)[0]
        self.assertAlmostEqual(
            test_result, 990.7401, delta=0.01,
            msg="Test result does not match expected value")

    def test_calc_ompc(self):
        """Test `_calc_ompc`.

        Use one set of inputs to test the estimation of total organic
        matter against a result calculated by hand.

        Raises:
            AssertionError if the raster created by `_calc_ompc`
                contains more than one unique value
            AssertionError if the value returned by `_calc_ompc' is
                not within 0.0001 of the value calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        som1c_2_path = os.path.join(self.workspace_dir, 'som1c_2.tif')
        som2c_2_path = os.path.join(self.workspace_dir, 'som2c_2.tif')
        som3c_path = os.path.join(self.workspace_dir, 'som3c.tif')
        bulk_d_path = os.path.join(self.workspace_dir, 'bulkd.tif')
        edepth_path = os.path.join(self.workspace_dir, 'edepth.tif')

        create_constant_raster(som1c_2_path, 42.109)
        create_constant_raster(som2c_2_path, 959.1091)
        create_constant_raster(som3c_path, 588.0574)
        create_constant_raster(bulk_d_path, 1.5)
        create_constant_raster(edepth_path, 0.2)

        ompc_path = os.path.join(self.workspace_dir, 'ompc.tif')

        forage._calc_ompc(
            som1c_2_path, som2c_2_path, som3c_path, bulk_d_path, edepth_path,
            ompc_path)

        # assert the value in the raster `ompc_path` is equal to value
        # calculated by hand
        result_set = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                ompc_path):
            result_set.update(numpy.unique(raster_block))
        self.assertEqual(
            len(result_set), 1,
            msg="One unique value expected in organic matter raster")
        test_result = list(result_set)[0]
        self.assertAlmostEqual(
            test_result, 0.913304, delta=0.0001,
            msg="Test result does not match expected value")

    def test_calc_afiel(self):
        """Test `_calc_afiel`.

        Use one set of inputs to test the estimation of field capacity against
        a result calculated by hand.

        Raises:
            AssertionError if the raster created by `_calc_afiel`
                contains more than one unique value
            AssertionError if the value returned by `_calc_afiel' is
                not within 0.0001 of the value calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        sand_path = os.path.join(self.workspace_dir, 'sand.tif')
        silt_path = os.path.join(self.workspace_dir, 'silt.tif')
        clay_path = os.path.join(self.workspace_dir, 'clay.tif')
        ompc_path = os.path.join(self.workspace_dir, 'ompc.tif')
        bulkd_path = os.path.join(self.workspace_dir, 'bulkd.tif')

        create_constant_raster(sand_path, 0.39)
        create_constant_raster(silt_path, 0.41)
        create_constant_raster(clay_path, 0.2)
        create_constant_raster(ompc_path, 0.913304)
        create_constant_raster(bulkd_path, 1.5)

        afiel_path = os.path.join(self.workspace_dir, 'afiel.tif')

        forage._calc_afiel(
            sand_path, silt_path, clay_path, ompc_path, bulkd_path, afiel_path)

        # assert the value in the raster `afiel_path` is equal to value
        # calculated by hand
        result_set = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                afiel_path):
            result_set.update(numpy.unique(raster_block))
        self.assertEqual(
            len(result_set), 1,
            msg="One unique value expected in field capacity raster")
        test_result = list(result_set)[0]
        self.assertAlmostEqual(
            test_result, 0.30895, delta=0.0001,
            msg="Test result does not match expected value")

    def test_calc_awilt(self):
        """Test `_calc_awilt`.

        Use one set of inputs to test the estimation of wilting point against
        a result calculated by hand.

        Raises:
            AssertionError if the raster created by `_calc_awilt`
                contains more than one unique value
            AssertionError if the value returned by `_calc_awilt' is
                not within 0.0001 of the value calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        sand_path = os.path.join(self.workspace_dir, 'sand.tif')
        silt_path = os.path.join(self.workspace_dir, 'silt.tif')
        clay_path = os.path.join(self.workspace_dir, 'clay.tif')
        ompc_path = os.path.join(self.workspace_dir, 'ompc.tif')
        bulkd_path = os.path.join(self.workspace_dir, 'bulkd.tif')

        create_constant_raster(sand_path, 0.39)
        create_constant_raster(silt_path, 0.41)
        create_constant_raster(clay_path, 0.2)
        create_constant_raster(ompc_path, 0.913304)
        create_constant_raster(bulkd_path, 1.5)

        awilt_path = os.path.join(self.workspace_dir, 'awilt.tif')

        forage._calc_awilt(
            sand_path, silt_path, clay_path, ompc_path, bulkd_path, awilt_path)

        # assert the value in the raster `awilt_path` is equal to value
        # calculated by hand
        result_set = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                awilt_path):
            result_set.update(numpy.unique(raster_block))
        self.assertEqual(
            len(result_set), 1,
            msg="One unique value expected in wilting point raster")
        test_result = list(result_set)[0]
        self.assertAlmostEqual(
            test_result, 0.201988, delta=0.0001,
            msg="Test result does not match expected value")

    def test_afiel_awilt(self):
        """Test `_afiel_awilt`.

        Use the function `_afiel_awilt` to calculate field capacity and wilting
        point from randomly generated inputs. Test that calculated field
        capacity and calculated wilting point are in the range [0.01, 0.9].
        Introduce nodata values into the input rasters and test that calculated
        field capacity and wilting point remain in the range [0.01, 0.9].
        Test the function with known values against results calculated by hand.

        Raises:
            AssertionError if a result from random inputs is outside the
                known possible range given the range of inputs
            AssertionError if a result from known inputs is outside the range
                [known result += 0.0001]

        Returns:
            None
        """
        from natcap.invest import forage

        site_param_table = {1: {'edepth': 0.2}}
        pp_reg = {
            'afiel_1_path': os.path.join(self.workspace_dir, 'afiel_1.tif'),
            'afiel_2_path': os.path.join(self.workspace_dir, 'afiel_2.tif'),
            'afiel_3_path': os.path.join(self.workspace_dir, 'afiel_3.tif'),
            'afiel_4_path': os.path.join(self.workspace_dir, 'afiel_4.tif'),
            'afiel_5_path': os.path.join(self.workspace_dir, 'afiel_5.tif'),
            'afiel_6_path': os.path.join(self.workspace_dir, 'afiel_6.tif'),
            'afiel_7_path': os.path.join(self.workspace_dir, 'afiel_7.tif'),
            'afiel_8_path': os.path.join(self.workspace_dir, 'afiel_8.tif'),
            'afiel_9_path': os.path.join(self.workspace_dir, 'afiel_9.tif'),
            'awilt_1_path': os.path.join(self.workspace_dir, 'awilt_1.tif'),
            'awilt_2_path': os.path.join(self.workspace_dir, 'awilt_2.tif'),
            'awilt_3_path': os.path.join(self.workspace_dir, 'awilt_3.tif'),
            'awilt_4_path': os.path.join(self.workspace_dir, 'awilt_4.tif'),
            'awilt_5_path': os.path.join(self.workspace_dir, 'awilt_5.tif'),
            'awilt_6_path': os.path.join(self.workspace_dir, 'awilt_6.tif'),
            'awilt_7_path': os.path.join(self.workspace_dir, 'awilt_7.tif'),
            'awilt_8_path': os.path.join(self.workspace_dir, 'awilt_8.tif'),
            'awilt_9_path': os.path.join(self.workspace_dir, 'awilt_9.tif')
            }

        site_index_path = os.path.join(self.workspace_dir, 'site_index.tif')
        som1c_2_path = os.path.join(self.workspace_dir, 'som1c_2.tif')
        som2c_2_path = os.path.join(self.workspace_dir, 'som2c_2.tif')
        som3c_path = os.path.join(self.workspace_dir, 'som3c.tif')
        sand_path = os.path.join(self.workspace_dir, 'sand.tif')
        silt_path = os.path.join(self.workspace_dir, 'silt.tif')
        clay_path = os.path.join(self.workspace_dir, 'clay.tif')
        bulk_d_path = os.path.join(self.workspace_dir, 'bulkd.tif')

        create_random_raster(site_index_path, 1, 1)
        create_random_raster(som1c_2_path, 35., 55.)
        create_random_raster(som2c_2_path, 500., 1500.)
        create_random_raster(som3c_path, 300., 600.)
        create_random_raster(sand_path, 0., 0.5)
        create_random_raster(silt_path, 0., 0.5)
        create_random_raster(bulk_d_path, 0.8, 1.8)
        create_complementary_raster(sand_path, silt_path, clay_path)

        minimum_acceptable_value = 0.01
        maximum_acceptable_value = 0.7
        nodata_value = _TARGET_NODATA

        forage._afiel_awilt(
            site_index_path, site_param_table, som1c_2_path,
            som2c_2_path, som3c_path, sand_path, silt_path, clay_path,
            bulk_d_path, pp_reg)

        for key, path in pp_reg.iteritems():
            self.assert_all_values_in_raster_within_range(
                path, minimum_acceptable_value,
                maximum_acceptable_value, nodata_value)

        for input_raster in [
                site_index_path, som1c_2_path, som2c_2_path, som3c_path,
                sand_path, silt_path, clay_path, bulk_d_path]:
            insert_nodata_values_into_raster(input_raster, _TARGET_NODATA)
            forage._afiel_awilt(
                site_index_path, site_param_table, som1c_2_path,
                som2c_2_path, som3c_path, sand_path, silt_path, clay_path,
                bulk_d_path, pp_reg)
            for key, path in pp_reg.iteritems():
                self.assert_all_values_in_raster_within_range(
                    path, minimum_acceptable_value,
                    maximum_acceptable_value, nodata_value)

        # known inputs
        site_param_table = {1: {'edepth': 0.111}}
        create_random_raster(som1c_2_path, 40, 40)
        create_random_raster(som2c_2_path, 744, 744)
        create_random_raster(som3c_path, 444, 444)
        create_random_raster(sand_path, 0.4, 0.4)
        create_random_raster(silt_path, 0.1, 0.1)
        create_random_raster(bulk_d_path, 0.81, 0.81)
        create_complementary_raster(sand_path, silt_path, clay_path)

        known_afiel_1 = 0.47285
        known_awilt_1 = 0.32424
        known_afiel_4 = 0.47085
        known_awilt_6 = 0.32132
        tolerance = 0.0001

        forage._afiel_awilt(
            site_index_path, site_param_table, som1c_2_path,
            som2c_2_path, som3c_path, sand_path, silt_path, clay_path,
            bulk_d_path, pp_reg)

        self.assert_all_values_in_raster_within_range(
            pp_reg['afiel_1_path'], known_afiel_1 - tolerance,
            known_afiel_1 + tolerance, nodata_value)
        self.assert_all_values_in_raster_within_range(
            pp_reg['awilt_1_path'], known_awilt_1 - tolerance,
            known_awilt_1 + tolerance, nodata_value)
        self.assert_all_values_in_raster_within_range(
            pp_reg['afiel_4_path'], known_afiel_4 - tolerance,
            known_afiel_4 + tolerance, nodata_value)
        self.assert_all_values_in_raster_within_range(
            pp_reg['awilt_6_path'], known_awilt_6 - tolerance,
            known_awilt_6 + tolerance, nodata_value)

    def test_persistent_params(self):
        """Test `persistent_params`.

        Use the function `persistent_params` to calculate wc, eftext, p1co2_2,
        fps1s3, and fps2s3 from randomly generated inputs. Test that each of
        the calculated quantities are within the range [0, 1].  Introduce
        nodata values into the inputs and test that calculated values
        remain inside the specified ranges. Test the function with known inputs
        against values calculated by hand.

        Raises:
            AssertionError if a result from random inputs is outside the
                known possible range given the range of inputs
            AssertionError if a result from known inputs is outside the range
                [known result += 0.0001]

        Returns:
            None
        """
        from natcap.invest import forage

        site_param_table = {
            1: {
                'peftxa': numpy.random.uniform(0.15, 0.35),
                'peftxb': numpy.random.uniform(0.65, 0.85),
                'p1co2a_2': numpy.random.uniform(0.1, 0.2),
                'p1co2b_2': numpy.random.uniform(0.58, 0.78),
                'ps1s3_1': numpy.random.uniform(0.58, 0.78),
                'ps1s3_2': numpy.random.uniform(0.02, 0.04),
                'ps2s3_1': numpy.random.uniform(0.58, 0.78),
                'ps2s3_2': numpy.random.uniform(0.001, 0.005),
                'omlech_1': numpy.random.uniform(0.01, 0.05),
                'omlech_2': numpy.random.uniform(0.06, 0.18)},
                }

        pp_reg = {
            'afiel_1_path': os.path.join(self.workspace_dir, 'afiel_1.tif'),
            'awilt_1_path': os.path.join(self.workspace_dir, 'awilt.tif'),
            'wc_path': os.path.join(self.workspace_dir, 'wc.tif'),
            'eftext_path': os.path.join(self.workspace_dir, 'eftext.tif'),
            'p1co2_2_path': os.path.join(self.workspace_dir, 'p1co2_2.tif'),
            'fps1s3_path': os.path.join(self.workspace_dir, 'fps1s3.tif'),
            'fps2s3_path': os.path.join(self.workspace_dir, 'fps2s3.tif'),
            'orglch_path': os.path.join(self.workspace_dir, 'orglch.tif'),
        }

        site_index_path = os.path.join(self.workspace_dir, 'site_index.tif')
        sand_path = os.path.join(self.workspace_dir, 'sand.tif')
        clay_path = os.path.join(self.workspace_dir, 'clay.tif')

        create_random_raster(site_index_path, 1, 1)
        create_random_raster(sand_path, 0., 0.5)
        create_random_raster(clay_path, 0., 0.5)
        create_random_raster(pp_reg['afiel_1_path'], 0.5, 0.9)
        create_random_raster(pp_reg['awilt_1_path'], 0.01, 0.49)

        acceptable_range_dict = {
            'wc_path': {
                'minimum_acceptable_value': 0.01,
                'maximum_acceptable_value': 0.89,
                'nodata_value': _TARGET_NODATA,
                },
            'eftext_path': {
                'minimum_acceptable_value': 0.15,
                'maximum_acceptable_value': 0.775,
                'nodata_value': _IC_NODATA,
                },
            'p1co2_2_path': {
                'minimum_acceptable_value': 0.1,
                'maximum_acceptable_value': 0.59,
                'nodata_value': _IC_NODATA,
                },
            'fps1s3_path': {
                'minimum_acceptable_value': 0.58,
                'maximum_acceptable_value': 0.8,
                'nodata_value': _IC_NODATA,
                },
            'fps2s3_path': {
                'minimum_acceptable_value': 0.58,
                'maximum_acceptable_value': 0.7825,
                'nodata_value': _IC_NODATA,
                },
            'orglch_path': {
                'minimum_acceptable_value': 0.01,
                'maximum_acceptable_value': 0.14,
                'nodata_value': _IC_NODATA,
                },
        }

        forage._persistent_params(
            site_index_path, site_param_table, sand_path, clay_path, pp_reg)

        for path, ranges in acceptable_range_dict.iteritems():
            self.assert_all_values_in_raster_within_range(
                pp_reg[path], ranges['minimum_acceptable_value'],
                ranges['maximum_acceptable_value'],
                ranges['nodata_value'])

        for input_raster in [
                site_index_path, sand_path, clay_path]:
            insert_nodata_values_into_raster(input_raster, _TARGET_NODATA)
            forage._persistent_params(
                site_index_path, site_param_table, sand_path, clay_path,
                pp_reg)

            for path, ranges in acceptable_range_dict.iteritems():
                self.assert_all_values_in_raster_within_range(
                    pp_reg[path], ranges['minimum_acceptable_value'],
                    ranges['maximum_acceptable_value'],
                    ranges['nodata_value'])

        # known inputs
        site_param_table[1]['peftxa'] = 0.2
        site_param_table[1]['peftxb'] = 0.7
        site_param_table[1]['p1co2a_2'] = 0.18
        site_param_table[1]['p1co2b_2'] = 0.65
        site_param_table[1]['ps1s3_1'] = 0.59
        site_param_table[1]['ps1s3_2'] = 0.03
        site_param_table[1]['ps2s3_1'] = 0.61
        site_param_table[1]['ps2s3_2'] = 0.0022
        site_param_table[1]['omlech_1'] = 0.022
        site_param_table[1]['omlech_2'] = 0.12

        create_random_raster(sand_path, 0.22, 0.22)
        create_random_raster(clay_path, 0.22, 0.22)
        create_random_raster(pp_reg['afiel_1_path'], 0.56, 0.56)
        create_random_raster(pp_reg['awilt_1_path'], 0.4, 0.4)

        known_value_dict = {
            'wc_path': {
                'value': 0.16,
                'nodata_value': _TARGET_NODATA,
                },
            'eftext_path': {
                'value': 0.354,
                'nodata_value': _IC_NODATA,
                },
            'p1co2_2_path': {
                'value': 0.323,
                'nodata_value': _IC_NODATA,
                },
            'fps1s3_path': {
                'value': 0.5966,
                'nodata_value': _IC_NODATA,
                },
            'fps2s3_path': {
                'value': 0.61048,
                'nodata_value': _IC_NODATA,
                },
            'orglch_path': {
                'value': 0.0484,
                'nodata_value': _IC_NODATA,
                },
        }
        tolerance = 0.0001

        forage._persistent_params(
            site_index_path, site_param_table, sand_path, clay_path,
            pp_reg)

        for path, values in known_value_dict.iteritems():
            self.assert_all_values_in_raster_within_range(
                pp_reg[path], values['value'] - tolerance,
                values['value'] + tolerance, values['nodata_value'])

    def test_aboveground_ratio(self):
        """Test `_aboveground_ratio`.

        Use the function `_aboveground_ratio` to calculate the C/N or P
        ratio of decomposing aboveground material from random inputs. Test
        that the calculated ratio, agdrat, is within the range [1, 150].
        Introduce nodata values into the inputs and test that calculated
        agdrat remains inside the range [1, 150]. Calculate aboveground
        ratio from known inputs and compare to result calculated by hand.

        Raises:
            AssertionError if agdrat from random inputs is outside the
                known possible range given the range of inputs
            AssertionError if agdrat from known inputs is outside the range
                [known result += 0.0001]

        Returns:
            None
        """
        from natcap.invest import forage

        array_shape = (10, 10)
        tolerance = 0.0001

        tca = numpy.random.uniform(300, 700, array_shape)
        anps = numpy.random.uniform(1, numpy.amin(tca), array_shape)
        pcemic_1 = numpy.random.uniform(12, 20, array_shape)
        pcemic_2 = numpy.random.uniform(3, 11, array_shape)
        pcemic_3 = numpy.random.uniform(0.001, 0.1, array_shape)

        minimum_acceptable_agdrat = 2.285
        maximum_acceptable_agdrat = numpy.amax(pcemic_1)
        agdrat_nodata = _TARGET_NODATA

        agdrat = forage._aboveground_ratio(
            anps, tca, pcemic_1, pcemic_2, pcemic_3)

        self.assert_all_values_in_array_within_range(
            agdrat, minimum_acceptable_agdrat, maximum_acceptable_agdrat,
            agdrat_nodata)

        for input_array in [anps, tca]:
            insert_nodata_values_into_array(input_array, _TARGET_NODATA)
            agdrat = forage._aboveground_ratio(
                anps, tca, pcemic_1, pcemic_2, pcemic_3)

            self.assert_all_values_in_array_within_range(
                agdrat, minimum_acceptable_agdrat, maximum_acceptable_agdrat,
                agdrat_nodata)
        for input_array in [pcemic_1, pcemic_2, pcemic_3]:
            insert_nodata_values_into_array(input_array, _IC_NODATA)
            agdrat = forage._aboveground_ratio(
                anps, tca, pcemic_1, pcemic_2, pcemic_3)

            self.assert_all_values_in_array_within_range(
                agdrat, minimum_acceptable_agdrat, maximum_acceptable_agdrat,
                agdrat_nodata)

        # known inputs: econt > pcemic_3
        tca = 413
        anps = 229
        pcemic_1 = 17.4
        pcemic_2 = 3.2
        pcemic_3 = 0.04

        known_agdrat = 3.2
        point_agdrat = agdrat_point(anps, tca, pcemic_1, pcemic_2, pcemic_3)
        self.assertAlmostEqual(known_agdrat, point_agdrat)

        tca_ar = numpy.full(array_shape, tca)
        anps_ar = numpy.full(array_shape, anps)
        pcemic_1_ar = numpy.full(array_shape, pcemic_1)
        pcemic_2_ar = numpy.full(array_shape, pcemic_2)
        pcemic_3_ar = numpy.full(array_shape, pcemic_3)

        agdrat = forage._aboveground_ratio(
            anps_ar, tca_ar, pcemic_1_ar, pcemic_2_ar, pcemic_3_ar)
        self.assert_all_values_in_array_within_range(
            agdrat, point_agdrat - tolerance, point_agdrat + tolerance,
            agdrat_nodata)

        # known inputs: econt < pcemic_3
        tca = 413.
        anps = 100.
        pcemic_1 = 17.4
        pcemic_2 = 3.2
        pcemic_3 = 0.11
        point_agdrat = agdrat_point(anps, tca, pcemic_1, pcemic_2, pcemic_3)

        tca_ar = numpy.full(array_shape, tca)
        anps_ar = numpy.full(array_shape, anps)
        pcemic_1_ar = numpy.full(array_shape, pcemic_1)
        pcemic_2_ar = numpy.full(array_shape, pcemic_2)
        pcemic_3_ar = numpy.full(array_shape, pcemic_3)

        agdrat = forage._aboveground_ratio(
            anps_ar, tca_ar, pcemic_1_ar, pcemic_2_ar, pcemic_3_ar)
        self.assert_all_values_in_array_within_range(
            agdrat, point_agdrat - tolerance, point_agdrat + tolerance,
            agdrat_nodata)

    def test_structural_ratios(self):
        """Test `_structural_ratios`.

        Use the function `_structural_ratios` to calculate rnewas_1_1,
        rnewas_1_2, rnewas_2_1, rnewas_2_2, rnewbs_1_2, and rnewbs_2_2 from
        randomly generated inputs. Test that each of the calculated quantities
        are within the range [1, 1500].  Introduce nodata values into the
        inputs and test that calculated values remain inside the specified
        ranges.

        Raises:
            AssertionError if rnewas_1_1 is outside the range [1, 1500]
            AssertionError if rnewas_1_2 is outside the range [1, 1500]
            AssertionError if rnewas_2_1 is outside the range [1, 1500]
            AssertionError if rnewas_2_2 is outside the range [1, 1500]
            AssertionError if rnewbs_1_1 is outside the range [1, 1500]
            AssertionError if rnewbs_2_2 is outside the range [1, 1500]
            AssertionError if rnewbs_2_1 is outside the range [1, 1500]
            AssertionError if rnewbs_2_2 is outside the range [1, 1500]

        Returns:
            None
        """
        from natcap.invest import forage

        site_param_table = {
            1: {
                'pcemic1_2_1': numpy.random.uniform(5, 12),
                'pcemic1_1_1': numpy.random.uniform(13, 23),
                'pcemic1_3_1': numpy.random.uniform(0.01, 0.05),
                'pcemic2_2_1': numpy.random.uniform(5, 12),
                'pcemic2_1_1': numpy.random.uniform(13, 23),
                'pcemic2_3_1': numpy.random.uniform(0.01, 0.05),
                'rad1p_1_1': numpy.random.uniform(8, 16),
                'rad1p_2_1': numpy.random.uniform(2, 5),
                'rad1p_3_1': numpy.random.uniform(2, 5),
                'varat1_1_1': numpy.random.uniform(12, 16),
                'varat22_1_1': numpy.random.uniform(15, 25),

                'pcemic1_2_2': numpy.random.uniform(90, 110),
                'pcemic1_1_2': numpy.random.uniform(170, 230),
                'pcemic1_3_2': numpy.random.uniform(0.0005, 0.0025),
                'pcemic2_2_2': numpy.random.uniform(75, 125),
                'pcemic2_1_2': numpy.random.uniform(200, 300),
                'pcemic2_3_2': numpy.random.uniform(0.0005, 0.0025),
                'rad1p_1_2': numpy.random.uniform(200, 300),
                'rad1p_2_2': numpy.random.uniform(3, 7),
                'rad1p_3_2': numpy.random.uniform(50, 150),
                'varat1_1_2': numpy.random.uniform(125, 175),
                'varat22_1_2': numpy.random.uniform(350, 450)},
                }

        sv_reg = {
            'strucc_1_path': os.path.join(self.workspace_dir, 'strucc_1.tif'),
            'struce_1_1_path': os.path.join(
                self.workspace_dir, 'struce_1_1.tif'),
            'struce_1_2_path': os.path.join(
                self.workspace_dir, 'struce_1_2.tif'),

        }
        site_index_path = os.path.join(self.workspace_dir, 'site_index.tif')
        create_random_raster(site_index_path, 1, 1)
        create_random_raster(sv_reg['strucc_1_path'], 120, 1800)
        create_random_raster(sv_reg['struce_1_1_path'], 0.5, 10)
        create_random_raster(sv_reg['struce_1_2_path'], 0.1, 0.50)

        pp_reg = {
            'rnewas_1_1_path': os.path.join(
                self.workspace_dir, 'rnewas_1_1.tif'),
            'rnewas_1_2_path': os.path.join(
                self.workspace_dir, 'rnewas_1_2.tif'),
            'rnewas_2_1_path': os.path.join(
                self.workspace_dir, 'rnewas_2_1.tif'),
            'rnewas_2_2_path': os.path.join(
                self.workspace_dir, 'rnewas_2_2.tif'),
            'rnewbs_1_1_path': os.path.join(
                self.workspace_dir, 'rnewbs_1_1.tif'),
            'rnewbs_1_2_path': os.path.join(
                self.workspace_dir, 'rnewbs_1_2.tif'),
            'rnewbs_2_1_path': os.path.join(
                self.workspace_dir, 'rnewbs_2_1.tif'),
            'rnewbs_2_2_path': os.path.join(
                self.workspace_dir, 'rnewbs_2_2.tif'),
        }

        minimum_acceptable_value = 1
        maximum_acceptable_value = 1500
        nodata_value = _TARGET_NODATA

        forage._structural_ratios(
            site_index_path, site_param_table, sv_reg, pp_reg)

        for key, path in pp_reg.iteritems():
            self.assert_all_values_in_raster_within_range(
                path, minimum_acceptable_value,
                maximum_acceptable_value, nodata_value)

        for input_raster in [
                site_index_path, sv_reg['strucc_1_path'],
                sv_reg['struce_1_1_path'], sv_reg['struce_1_2_path']]:
            insert_nodata_values_into_raster(input_raster, _TARGET_NODATA)
            forage._structural_ratios(
                site_index_path, site_param_table, sv_reg, pp_reg)

            for key, path in pp_reg.iteritems():
                self.assert_all_values_in_raster_within_range(
                    path, minimum_acceptable_value,
                    maximum_acceptable_value, nodata_value)

    def test_yearly_tasks(self):
        """Test `_yearly_tasks`.

        Use `_yearly_tasks` to calculate annual precipitation and annual
        atmospheric N deposition from random inputs. Test that the function
        fails if fewer than 12 months of precipitation are supplied. Test that
        the function fails if the dates of precipitation inputs do not fill
        the 12 months surrounding the current month. Test that annual
        precipitation falls inside the range [0, 72]. Test that annual
        atmospheric N deposition falls inside the range [0, 37]. Introduce
        nodata values into input rasters and test that atmospheric N deposition
        remains within the specified range.

        Raises:
            AssertionError if `_yearly_tasks` does not fail with fewer than 12
                months of precipitation rasters supplied
            AssertionError if `_yearly_tasks` does not fail with precipitation
                rasters supplied not within 12 months of current month
            AssertionError if a result from random inputs is outside the
                known possible range given the range of inputs

        Returns:
            None
        """
        from natcap.invest import forage

        month_index = numpy.random.randint(0, 100)
        site_index_path = os.path.join(self.workspace_dir, 'site_index.tif')
        site_param_table = {
            1: {
                'epnfa_1': numpy.random.uniform(0, 1),
                'epnfa_2': numpy.random.uniform(0, 0.5),
                }
            }

        complete_aligned_inputs = {
            'precip_{}'.format(month): os.path.join(
                self.workspace_dir, 'precip_{}.tif'.format(month)) for
            month in xrange(month_index, month_index + 12)
        }

        year_reg = {
            'annual_precip_path': os.path.join(
                self.workspace_dir, 'annual_precip.tif'),
            'baseNdep_path': os.path.join(self.workspace_dir, 'baseNdep.tif')
        }

        create_random_raster(site_index_path, 1, 1)
        for key, path in complete_aligned_inputs.iteritems():
            create_random_raster(path, 0, 6)

        # fewer than 12 months of precip rasters
        modified_inputs = complete_aligned_inputs.copy()
        removed_key = modified_inputs.pop('precip_{}'.format(
            numpy.random.randint(month_index, month_index + 12)))
        with self.assertRaises(ValueError):
            forage._yearly_tasks(
                site_index_path, site_param_table, modified_inputs,
                month_index, year_reg)

        # 12 months of precip rasters supplied, but outside 12 month window of
        # current month
        modified_inputs['precip_{}'.format(month_index + 13)] = os.path.join(
            'precip_{}.tif'.format(month_index + 13))
        with self.assertRaises(ValueError):
            forage._yearly_tasks(
                site_index_path, site_param_table, modified_inputs,
                month_index, year_reg)

        # complete intact inputs
        minimum_acceptable_annual_precip = 0
        maximum_acceptabe_annual_precip = 72
        precip_nodata = _TARGET_NODATA

        minimum_acceptable_Ndep = 0
        maximum_acceptable_Ndep = 37
        Ndep_nodata = _TARGET_NODATA

        forage._yearly_tasks(
            site_index_path, site_param_table, complete_aligned_inputs,
            month_index, year_reg)
        self.assert_all_values_in_raster_within_range(
            year_reg['annual_precip_path'], minimum_acceptable_annual_precip,
            maximum_acceptabe_annual_precip, precip_nodata)
        self.assert_all_values_in_raster_within_range(
            year_reg['baseNdep_path'], minimum_acceptable_Ndep,
            maximum_acceptable_Ndep, Ndep_nodata)

        input_raster_list = [site_index_path] + [
            path for key, path in complete_aligned_inputs.iteritems()]
        for input_raster in input_raster_list:
            insert_nodata_values_into_raster(input_raster, _TARGET_NODATA)
            forage._yearly_tasks(
                site_index_path, site_param_table, complete_aligned_inputs,
                month_index, year_reg)
            self.assert_all_values_in_raster_within_range(
                year_reg['annual_precip_path'],
                minimum_acceptable_annual_precip,
                maximum_acceptabe_annual_precip, precip_nodata)
            self.assert_all_values_in_raster_within_range(
                year_reg['baseNdep_path'], minimum_acceptable_Ndep,
                maximum_acceptable_Ndep, Ndep_nodata)

    def test_reference_evapotranspiration(self):
        """Test `_reference_evapotranspiration`.

        Use the function `_reference_evapotranspiration` to calculate reference
        evapotranspiration (ET) from random inputs. Test that the calculated
        reference ET is within the range [0, 32]. Introduce nodata values into
        the inputs and test that the result remains inside the range [0, 31].
        Test the function with known inputs against a value calculated by hand.

        Raises:
            AssertionError if evapotranspiration from random inputs is outside
                the known possible range given the range of inputs
            AssertionError if evapotranspiration from known inputs is outside
                the range [known result += 0.0001]

        Returns:
            None
        """
        from natcap.invest import forage

        max_temp_path = os.path.join(self.workspace_dir, 'max_temp.tif')
        min_temp_path = os.path.join(self.workspace_dir, 'min_temp.tif')
        shwave_path = os.path.join(self.workspace_dir, 'shwave.tif')
        fwloss_4_path = os.path.join(self.workspace_dir, 'fwloss_4.tif')

        create_random_raster(max_temp_path, 21, 40)
        create_random_raster(min_temp_path, -20, 20)
        create_random_raster(shwave_path, 0, 1125)
        create_random_raster(fwloss_4_path, 0, 1)

        pevap_path = os.path.join(self.workspace_dir, 'pevap.tif')

        minimum_acceptable_ET = 0
        maximum_acceptable_ET = 32
        ET_nodata = _TARGET_NODATA

        forage._reference_evapotranspiration(
            max_temp_path, min_temp_path, shwave_path, fwloss_4_path,
            pevap_path)

        self.assert_all_values_in_raster_within_range(
            pevap_path, minimum_acceptable_ET, maximum_acceptable_ET,
            ET_nodata)

        insert_nodata_values_into_raster(max_temp_path, _IC_NODATA)
        insert_nodata_values_into_raster(min_temp_path, _IC_NODATA)
        insert_nodata_values_into_raster(shwave_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(fwloss_4_path, _IC_NODATA)

        forage._reference_evapotranspiration(
            max_temp_path, min_temp_path, shwave_path, fwloss_4_path,
            pevap_path)

        self.assert_all_values_in_raster_within_range(
            pevap_path, minimum_acceptable_ET, maximum_acceptable_ET,
            ET_nodata)

        # known inputs
        create_random_raster(max_temp_path, 23, 23)
        create_random_raster(min_temp_path, -2, -2)
        create_random_raster(shwave_path, 880, 880)
        create_random_raster(fwloss_4_path, 0.6, 0.6)

        known_ET = 9.5465
        tolerance = 0.0001

        forage._reference_evapotranspiration(
            max_temp_path, min_temp_path, shwave_path, fwloss_4_path,
            pevap_path)

        self.assert_all_values_in_raster_within_range(
            pevap_path, known_ET - tolerance, known_ET + tolerance,
            ET_nodata)

    def test_potential_production(self):
        """Test `_potential_production`.

        Use the function `_potential_production` to calculate h2ogef_1 and
        total potential production from random inputs. Test that h2ogef_1, the
        limiting factor of water availability on growth, is inside the range
        [0.009, 1]. Test that total production is inside the range [0, 675].
        Introduce nodata values into inputs and test that h2ogef_1 and
        potential production remain inside the specified ranges.

        Raises:
            AssertionError if h2ogef_1 is outside the range [0.009, 1]
            AssertionError if total potential production is outside the range
                [0, 675]

        Returns:
            None
        """
        from natcap.invest import forage

        month_index = 10
        current_month = 6
        pft_id_set = set([1, 2])

        aligned_inputs = {
            'site_index': os.path.join(self.workspace_dir, 'site_index.tif'),
            'max_temp_{}'.format(current_month): os.path.join(
                self.workspace_dir, 'max_temp.tif'),
            'min_temp_{}'.format(current_month): os.path.join(
                self.workspace_dir, 'min_temp.tif'),
            'precip_{}'.format(month_index): os.path.join(
                self.workspace_dir, 'precip.tif'),
        }
        create_random_raster(aligned_inputs['site_index'], 1, 1)
        create_random_raster(
            aligned_inputs['max_temp_{}'.format(current_month)], 10, 30)
        create_random_raster(
            aligned_inputs['min_temp_{}'.format(current_month)], -10, 9)
        create_random_raster(
            aligned_inputs['precip_{}'.format(month_index)], 0, 6)

        for pft_i in pft_id_set:
            aligned_inputs['pft_{}'.format(pft_i)] = os.path.join(
                self.workspace_dir, 'pft_{}.tif'.format(pft_i))
            create_random_raster(
                aligned_inputs['pft_{}'.format(pft_i)], 0, 1)

        site_param_table = {
            1: {
                'pmxbio': numpy.random.uniform(500, 700),
                'pmxtmp': numpy.random.uniform(-0.0025, 0),
                'pmntmp': numpy.random.uniform(0, 0.01),
                'fwloss_4': numpy.random.uniform(0, 1),
                'pprpts_1': numpy.random.uniform(0, 1),
                'pprpts_2': numpy.random.uniform(0.5, 1.5),
                'pprpts_3': numpy.random.uniform(0, 1),
            }
        }
        veg_trait_table = {}
        for pft_i in pft_id_set:
            veg_trait_table[pft_i] = {
                'ppdf_1': numpy.random.uniform(10, 30),
                'ppdf_2': numpy.random.uniform(31, 50),
                'ppdf_3': numpy.random.uniform(0, 1),
                'ppdf_4': numpy.random.uniform(0, 10),
                'biok5': numpy.random.uniform(0, 2000),
                'prdx_1': numpy.random.uniform(0.1, 0.6),
                'growth_months': ['3', '4', '5', '6'],
            }

        sv_reg = {
            'strucc_1_path': os.path.join(self.workspace_dir, 'strucc_1.tif'),
        }
        create_random_raster(sv_reg['strucc_1_path'], 0, 200)
        for pft_i in pft_id_set:
            sv_reg['aglivc_{}_path'.format(pft_i)] = os.path.join(
                self.workspace_dir, 'aglivc_{}.tif'.format(pft_i))
            create_random_raster(sv_reg['aglivc_{}_path'.format(pft_i)], 0, 50)
            sv_reg['stdedc_{}_path'.format(pft_i)] = os.path.join(
                self.workspace_dir, 'stdedc_{}.tif'.format(pft_i))
            create_random_raster(sv_reg['stdedc_{}_path'.format(pft_i)], 0, 50)
            sv_reg['avh2o_1_{}_path'.format(pft_i)] = os.path.join(
                self.workspace_dir, 'avh2o_1_{}.tif'.format(pft_i))
            create_random_raster(
                sv_reg['avh2o_1_{}_path'.format(pft_i)], 0, 3.5)

        pp_reg = {
            'wc_path': os.path.join(self.workspace_dir, 'wc.tif')
        }
        create_random_raster(pp_reg['wc_path'], 0.01, 0.9)

        month_reg = {}
        for pft_i in pft_id_set:
            month_reg['h2ogef_1_{}'.format(pft_i)] = os.path.join(
                self.workspace_dir, 'h2ogef_1_{}.tif'.format(pft_i))
            month_reg['tgprod_pot_prod_{}'.format(pft_i)] = os.path.join(
                self.workspace_dir, 'tgprod_pot_prod_{}.tif'.format(pft_i))

        minimum_acceptable_h2ogef_1 = 0.009
        maximum_acceptable_h2ogef_1 = 1

        minimum_acceptable_potential_production = 0
        maximum_acceptable_potential_production = 675

        forage._potential_production(
            aligned_inputs, site_param_table, current_month, month_index,
            pft_id_set, veg_trait_table, sv_reg, pp_reg, month_reg)

        for pft_i in pft_id_set:
            self.assert_all_values_in_raster_within_range(
                month_reg['h2ogef_1_{}'.format(pft_i)],
                minimum_acceptable_h2ogef_1,
                maximum_acceptable_h2ogef_1, _TARGET_NODATA)
            self.assert_all_values_in_raster_within_range(
                month_reg['tgprod_pot_prod_{}'.format(pft_i)],
                minimum_acceptable_potential_production,
                maximum_acceptable_potential_production, _TARGET_NODATA)

        insert_nodata_values_into_raster(
            aligned_inputs['site_index'], _TARGET_NODATA)
        insert_nodata_values_into_raster(
            aligned_inputs['max_temp_{}'.format(current_month)],
            -9999)
        insert_nodata_values_into_raster(
            aligned_inputs['min_temp_{}'.format(current_month)],
            _IC_NODATA)
        insert_nodata_values_into_raster(
            aligned_inputs['pft_{}'.format(list(pft_id_set)[0])],
            _TARGET_NODATA)
        insert_nodata_values_into_raster(
            sv_reg['aglivc_{}_path'.format(list(pft_id_set)[0])], _SV_NODATA)
        insert_nodata_values_into_raster(
            sv_reg['avh2o_1_{}_path'.format(list(pft_id_set)[1])],
            _TARGET_NODATA)
        insert_nodata_values_into_raster(pp_reg['wc_path'], _TARGET_NODATA)

        forage._potential_production(
            aligned_inputs, site_param_table, current_month, month_index,
            pft_id_set, veg_trait_table, sv_reg, pp_reg, month_reg)

        for pft_i in pft_id_set:
            self.assert_all_values_in_raster_within_range(
                month_reg['h2ogef_1_{}'.format(pft_i)],
                minimum_acceptable_h2ogef_1,
                maximum_acceptable_h2ogef_1, _TARGET_NODATA)
            self.assert_all_values_in_raster_within_range(
                month_reg['tgprod_pot_prod_{}'.format(pft_i)],
                minimum_acceptable_potential_production,
                maximum_acceptable_potential_production, _TARGET_NODATA)

    def test_calc_favail_P(self):
        """Test `_calc_favail_P`.

        Use the function `_calc_favail_P` to calculate the intermediate
        parameter favail_P from random inputs.  Test that favail_P is
        inside the range [0, 1]. Introduce nodata values into inputs and test
        that favail_P remains inside the range [0, 1]. Test the function with
        known inputs against values calculated by hand.

        Raises:
            AssertionError if favail_P from random inputs is outside the
                known possible range given the range of inputs
            AssertionError if favail_ from known inputs is outside the range
                [known result += 0.0001]

        Returns:
            None
        """
        from natcap.invest import forage

        sv_reg = {
            'minerl_1_1_path': os.path.join(
                self.workspace_dir, 'minerl_1_1.tif')
        }
        param_val_dict = {
            'favail_4': os.path.join(self.workspace_dir, 'favail_4.tif'),
            'favail_5': os.path.join(self.workspace_dir, 'favail_5.tif'),
            'favail_6': os.path.join(self.workspace_dir, 'favail_6.tif'),
            'favail_2': os.path.join(self.workspace_dir, 'favail_2.tif'),
        }

        create_random_raster(sv_reg['minerl_1_1_path'], 3, 8)
        create_random_raster(param_val_dict['favail_4'], 0, 1)
        create_random_raster(param_val_dict['favail_5'], 0, 1)
        create_random_raster(param_val_dict['favail_6'], 1, 3)

        forage._calc_favail_P(sv_reg, param_val_dict)

        minimum_acceptable_favail_P = 0
        maximum_acceptable_favail_P = 1

        self.assert_all_values_in_raster_within_range(
            param_val_dict['favail_2'],
            minimum_acceptable_favail_P,
            maximum_acceptable_favail_P, _IC_NODATA)

        insert_nodata_values_into_raster(sv_reg['minerl_1_1_path'], _SV_NODATA)
        for input_raster in [
                param_val_dict['favail_4'], param_val_dict['favail_5'],
                param_val_dict['favail_6']]:
            insert_nodata_values_into_raster(input_raster, _IC_NODATA)
            forage._calc_favail_P(sv_reg, param_val_dict)
            self.assert_all_values_in_raster_within_range(
                param_val_dict['favail_2'],
                minimum_acceptable_favail_P,
                maximum_acceptable_favail_P, _IC_NODATA)

        # known inputs
        create_random_raster(sv_reg['minerl_1_1_path'], 4.5, 4.5)
        create_random_raster(param_val_dict['favail_4'], 0.2, 0.2)
        create_random_raster(param_val_dict['favail_5'], 0.5, 0.5)
        create_random_raster(param_val_dict['favail_6'], 2.3, 2.3)

        known_favail_2 = 0.5
        tolerance = 0.0001

        forage._calc_favail_P(sv_reg, param_val_dict)
        self.assert_all_values_in_raster_within_range(
            param_val_dict['favail_2'],
            known_favail_2 - tolerance, known_favail_2 + tolerance,
            _IC_NODATA)

    def test_raster_list_sum(self):
        """Test `raster_list_sum`.

        Use the function `raster_list_sum` to calculate the sum across pixels
        in three rasters containing nodata.  Test that when
        nodata_remove=False, the result also contains nodata values. Test
        that when nodata_remove=True, nodata pixels are treated as zero.

        Raises:
            AssertionError if result raster does not contain nodata values
                in same position as input rasters

        Returns:
            None
        """
        from natcap.invest import forage

        num_rasters = numpy.random.randint(1, 10)
        raster_list = [
            os.path.join(self.workspace_dir, '{}.tif'.format(r)) for r in
            xrange(num_rasters)]

        for input_raster in raster_list:
            create_random_raster(input_raster, 1, 1)

        input_nodata = -999
        target_path = os.path.join(self.workspace_dir, 'result.tif')
        target_nodata = -9.99

        # input rasters include no nodata values
        forage.raster_list_sum(
            raster_list, input_nodata, target_path, target_nodata,
            nodata_remove=False)
        self.assert_all_values_in_raster_within_range(
            target_path, num_rasters, num_rasters, target_nodata)

        forage.raster_list_sum(
            raster_list, input_nodata, target_path, target_nodata,
            nodata_remove=True)
        self.assert_all_values_in_raster_within_range(
            target_path, num_rasters, num_rasters, target_nodata)

        # one input raster includes nodata values
        insert_nodata_values_into_raster(raster_list[0], input_nodata)

        forage.raster_list_sum(
            raster_list, input_nodata, target_path, target_nodata,
            nodata_remove=False)
        self.assert_all_values_in_raster_within_range(
            target_path, num_rasters, num_rasters, target_nodata)

        # assert that raster_list[0] and target_path include nodata
        # values in same locations
        input_including_nodata = gdal.OpenEx(raster_list[0])
        result_including_nodata = gdal.OpenEx(target_path)
        input_band = input_including_nodata.GetRasterBand(1)
        result_band = result_including_nodata.GetRasterBand(1)
        input_array = input_band.ReadAsArray()
        result_array = result_band.ReadAsArray()

        sum_input_mask = numpy.sum(input_array[input_array == input_nodata])
        sum_result_mask = numpy.sum(input_array[result_array == target_nodata])

        self.assertEqual(
            sum_input_mask, sum_result_mask,
            msg="Result raster must contain nodata values in same " +
            "position as input")

        input_band = None
        result_band = None
        input_including_nodata = None
        result_including_nodata = None

        forage.raster_list_sum(
            raster_list, input_nodata, target_path, target_nodata,
            nodata_remove=True)

        # assert that minimum value in target_path is num_rasters - 1
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                target_path):
            if len(raster_block[raster_block != target_nodata]) == 0:
                continue
            min_val = numpy.amin(
                raster_block[raster_block != target_nodata])
            self.assertGreaterEqual(
                min_val, (num_rasters - 1),
                msg="Raster appears to contain nodata values")

    def test_weighted_state_variable_sum(self):
        """Test `weighted_state_variable_sum`.

        Use the function `weighted_state_variable_sum` to calculate the
        weighted sum of a state variable across plant functional types. Test
        that the calculated sum matches values calculated by hand.

        Raises:
            AssertionError if the result calculated by
                `weighted_state_variable_sum` is outside the range
                [known result += 0.0001]

        Returns:
            None
        """
        from natcap.invest import forage

        sv = 'state_variable'
        pft_id_set = [2, 5, 7]
        percent_cover_dict = {
            pft_id_set[0]: 0.3,
            pft_id_set[1]: 0.001,
            pft_id_set[2]: 0.58,
        }
        sv_value_dict = {
            pft_id_set[0]: 20.,
            pft_id_set[1]: 300.84,
            pft_id_set[2]: 102.,
        }
        sv_reg = {}
        aligned_inputs = {}
        for pft_i in pft_id_set:
            aligned_inputs['pft_{}'.format(pft_i)] = os.path.join(
                self.workspace_dir, 'pft_{}.tif'.format(pft_i))
            create_constant_raster(
                aligned_inputs['pft_{}'.format(pft_i)],
                percent_cover_dict[pft_i])
            sv_reg['{}_{}_path'.format(sv, pft_i)] = os.path.join(
                self.workspace_dir, '{}_{}.tif'.format(sv, pft_i))
            create_constant_raster(
                sv_reg['{}_{}_path'.format(sv, pft_i)],
                sv_value_dict[pft_i])
        weighted_sum_path = os.path.join(
            self.workspace_dir, 'weighted_sum.tif')

        tolerance = 0.0001

        # known inputs
        known_weighted_sum = 65.46084
        forage.weighted_state_variable_sum(
            sv, sv_reg, aligned_inputs, pft_id_set, weighted_sum_path)
        self.assert_all_values_in_raster_within_range(
            weighted_sum_path, known_weighted_sum - tolerance,
            known_weighted_sum + tolerance, _TARGET_NODATA)

        # one pft has zero percent cover
        percent_cover_dict[pft_id_set[0]] = 0.
        create_constant_raster(
            aligned_inputs['pft_{}'.format(pft_id_set[0])],
            percent_cover_dict[pft_id_set[0]])

        known_weighted_sum = 59.46084
        forage.weighted_state_variable_sum(
            sv, sv_reg, aligned_inputs, pft_id_set, weighted_sum_path)
        self.assert_all_values_in_raster_within_range(
            weighted_sum_path, known_weighted_sum - tolerance,
            known_weighted_sum + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_raster(
            aligned_inputs['pft_{}'.format(pft_id_set[0])], _TARGET_NODATA)
        insert_nodata_values_into_raster(
            sv_reg['{}_{}_path'.format(sv, pft_id_set[2])], _SV_NODATA)

        forage.weighted_state_variable_sum(
            sv, sv_reg, aligned_inputs, pft_id_set, weighted_sum_path)

    def test_calc_available_nutrient(self):
        """Test `_calc_available_nutrient`.

        Use the function `_calc_available_nutrient` to calculate available
        nutrient from random results. Test that the calculated nutrient
        available is inside the range [0, 323]. Introduce nodata values into
        inputs and test that available nutrient remains inside the range
        [0, 323]. Test the function with known values against values calculated
        by hand.

        Raises:
            AssertionError if a result from random inputs is outside the
                known possible range given the range of inputs
            AssertionError if a result from known inputs is outside the range
                [known result += 0.0001]

        Returns:
            None
        """
        from natcap.invest import forage

        pft_i = numpy.random.randint(0, 4)
        pft_param_dict = {
            'snfxmx_1': numpy.random.uniform(0, 1),
            'nlaypg': numpy.random.randint(1, 10),
        }
        sv_reg = {
            'bglivc_{}_path'.format(pft_i): os.path.join(
                self.workspace_dir, 'bglivc_path.tif'),
            'crpstg_1_{}_path'.format(pft_i): os.path.join(
                self.workspace_dir, 'crpstg_1_{}.tif'.format(pft_i)),
            'crpstg_2_{}_path'.format(pft_i): os.path.join(
                self.workspace_dir, 'crpstg_2_{}.tif'.format(pft_i)),
        }
        for iel in [1, 2]:
            for lyr in xrange(1, 11):
                sv_reg['minerl_{}_{}_path'.format(lyr, iel)] = os.path.join(
                    self.workspace_dir, 'minerl_{}_{}.tif'.format(lyr, iel))
                create_random_raster(
                    sv_reg['minerl_{}_{}_path'.format(lyr, iel)], 0, 5)

        site_param_table = {
            1: {
                'rictrl': numpy.random.uniform(0.005, 0.02),
                'riint': numpy.random.uniform(0.6, 1),
                }
            }
        site_index_path = os.path.join(self.workspace_dir, 'site_index.tif')
        favail_path = os.path.join(self.workspace_dir, 'favail.tif')
        tgprod_path = os.path.join(self.workspace_dir, 'tgprod.tif')

        create_random_raster(site_index_path, 1, 1)
        create_random_raster(sv_reg['bglivc_{}_path'.format(pft_i)], 90, 180)
        create_random_raster(sv_reg['crpstg_1_{}_path'.format(pft_i)], 0, 3)
        create_random_raster(sv_reg['crpstg_2_{}_path'.format(pft_i)], 0, 1)
        create_random_raster(favail_path, 0, 1)
        create_random_raster(tgprod_path, 0, 675)

        eavail_path = os.path.join(self.workspace_dir, 'eavail.tif')

        minimum_acceptable_eavail = 0
        maximum_acceptable_evail = 323

        for iel in [1, 2]:
            forage._calc_available_nutrient(
                pft_i, iel, pft_param_dict, sv_reg, site_param_table,
                site_index_path, favail_path, tgprod_path, eavail_path)

            self.assert_all_values_in_raster_within_range(
                eavail_path, minimum_acceptable_eavail,
                maximum_acceptable_evail, _TARGET_NODATA)

        insert_nodata_values_into_raster(site_index_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(
            sv_reg['minerl_1_1_path'], _TARGET_NODATA)
        insert_nodata_values_into_raster(
            sv_reg['bglivc_{}_path'.format(pft_i)], _TARGET_NODATA)
        insert_nodata_values_into_raster(tgprod_path, _TARGET_NODATA)

        for iel in [1, 2]:
            forage._calc_available_nutrient(
                pft_i, iel, pft_param_dict, sv_reg, site_param_table,
                site_index_path, favail_path, tgprod_path, eavail_path)

            self.assert_all_values_in_raster_within_range(
                eavail_path, minimum_acceptable_eavail,
                maximum_acceptable_evail, _TARGET_NODATA)

        # known inputs
        create_random_raster(sv_reg['bglivc_{}_path'.format(pft_i)], 100, 100)
        create_random_raster(
            sv_reg['crpstg_1_{}_path'.format(pft_i)], 0.8, 0.8)
        create_random_raster(
            sv_reg['crpstg_2_{}_path'.format(pft_i)], 0.8, 0.8)
        create_random_raster(favail_path, 0.3, 0.3)
        create_random_raster(tgprod_path, 300, 300)

        for iel in [1, 2]:
            for lyr in xrange(1, 11):
                create_random_raster(
                    sv_reg['minerl_{}_{}_path'.format(lyr, iel)], 1, 1)

        pft_param_dict['snfxmx_1'] = 0.4
        pft_param_dict['nlaypg'] = 4

        site_param_table[1]['rictrl'] = 0.013
        site_param_table[1]['riint'] = 0.65

        known_N_avail = 49.9697
        known_P_avail = 1.9697
        tolerance = 0.0001

        iel = 1
        eavail_path = os.path.join(self.workspace_dir, 'eavail_N.tif')
        forage._calc_available_nutrient(
            pft_i, iel, pft_param_dict, sv_reg, site_param_table,
            site_index_path, favail_path, tgprod_path, eavail_path)

        self.assert_all_values_in_raster_within_range(
            eavail_path, known_N_avail - tolerance,
            known_N_avail + tolerance, _TARGET_NODATA)

        iel = 2
        eavail_path = os.path.join(self.workspace_dir, 'eavail_P.tif')
        forage._calc_available_nutrient(
            pft_i, iel, pft_param_dict, sv_reg, site_param_table,
            site_index_path, favail_path, tgprod_path, eavail_path)

        self.assert_all_values_in_raster_within_range(
            eavail_path, known_P_avail - tolerance,
            known_P_avail + tolerance, _TARGET_NODATA)

    def test_calc_nutrient_demand(self):
        """Test `_calc_nutrient_demand`.

        Use the function `_calc_nutrient_demand` to calculate demand for
        one nutrient by one plant functional type. Test that the calculated
        demand is in the range [???].  Introduce nodata values into inputs
        and test that calculated demand remains in the range [???]. Test
        against a result calculated by hand for known inputs.

        Raises:
            AssertionError if demand from random inputs is outside the
                known possible range given the range of inputs
            AssertionError if demand from known inputs is outside the range
                [known result += 0.0001]

        Returns:
            None
        """
        from natcap.invest import forage

        biomass_production_path = os.path.join(
            self.workspace_dir, 'biomass_production.tif')
        fraction_allocated_to_roots_path = os.path.join(
            self.workspace_dir, 'fraction_allocated_to_roots.tif')
        cercrp_min_above_path = os.path.join(
            self.workspace_dir, 'cercrp_min_above.tif')
        cercrp_min_below_path = os.path.join(
            self.workspace_dir, 'cercrp_min_below.tif')
        demand_path = os.path.join(
            self.workspace_dir, 'demand.tif')

        # run with random inputs
        create_random_raster(biomass_production_path, 0, 675)
        create_random_raster(fraction_allocated_to_roots_path, 0.01, 0.99)
        create_random_raster(cercrp_min_above_path, 8, 16)
        create_random_raster(cercrp_min_below_path, 8, 16)

        minimum_acceptable_demand = 0
        maximum_acceptable_demand = 33.75

        forage._calc_nutrient_demand(
            biomass_production_path, fraction_allocated_to_roots_path,
            cercrp_min_above_path, cercrp_min_below_path, demand_path)

        self.assert_all_values_in_raster_within_range(
            demand_path, minimum_acceptable_demand,
            maximum_acceptable_demand, _TARGET_NODATA)

        # insert nodata values into inputs
        insert_nodata_values_into_raster(
            biomass_production_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(
            fraction_allocated_to_roots_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(cercrp_min_above_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(cercrp_min_below_path, _TARGET_NODATA)

        forage._calc_nutrient_demand(
            biomass_production_path, fraction_allocated_to_roots_path,
            cercrp_min_above_path, cercrp_min_below_path, demand_path)

        self.assert_all_values_in_raster_within_range(
            demand_path, minimum_acceptable_demand,
            maximum_acceptable_demand, _TARGET_NODATA)

        # run with known inputs
        create_random_raster(biomass_production_path, 300, 300)
        create_random_raster(fraction_allocated_to_roots_path, 0.4, 0.4)
        create_random_raster(cercrp_min_above_path, 15, 15)
        create_random_raster(cercrp_min_below_path, 9, 9)

        known_demand = 10.1333
        tolerance = 0.0001

        forage._calc_nutrient_demand(
            biomass_production_path, fraction_allocated_to_roots_path,
            cercrp_min_above_path, cercrp_min_below_path, demand_path)

        self.assert_all_values_in_raster_within_range(
            demand_path, known_demand - tolerance, known_demand + tolerance,
            _TARGET_NODATA)

    def test_calc_provisional_fracrc(self):
        """Test `calc_provisional_fracrc`.

        Use the function `calc_provisional_fracrc` to calculate fracrc_p, the
        fraction of carbon allocated to roots. Test that fracrc_p calculated
        from random inputs is inside the valid range given the range
        of inputs. Introduce nodata values into inputs and test that fracrc_p
        remains inside the valid range. Test the function with known inputs
        against values calculated by hand.

        Raises:
            AssertionError if fracrc_p from random inputs is outside the range
                of valid values given the range of inputs
            AssertionError if fracrc_p from known inputs is not within 0.0001
                of the value calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        array_shape = (10, 10)

        annual_precip = numpy.random.uniform(22, 100, array_shape)
        frtcindx = numpy.random.randint(0, 2, array_shape)
        bgppa = numpy.random.uniform(100, 200, array_shape)
        bgppb = numpy.random.uniform(2, 12, array_shape)
        agppa = numpy.random.uniform(-40, -10, array_shape)
        agppb = numpy.random.uniform(2, 12, array_shape)
        cfrtcw_1 = numpy.random.uniform(0.4, 0.8, array_shape)
        cfrtcw_2 = numpy.random.uniform(0.01, 0.38, array_shape)
        cfrtcn_1 = numpy.random.uniform(0.4, 0.8, array_shape)
        cfrtcn_2 = numpy.random.uniform(0.01, 0.38, array_shape)

        minimum_acceptable_fracrc_p = 0.205
        maximum_acceptable_fracrc_p = 0.97297

        fracrc_p = forage.calc_provisional_fracrc(
            annual_precip, frtcindx, bgppa, bgppb, agppa, agppb,
            cfrtcw_1, cfrtcw_2, cfrtcn_1, cfrtcn_2)
        self.assert_all_values_in_array_within_range(
            fracrc_p, minimum_acceptable_fracrc_p,
            maximum_acceptable_fracrc_p, _TARGET_NODATA)

        insert_nodata_values_into_array(annual_precip, _TARGET_NODATA)
        fracrc_p = forage.calc_provisional_fracrc(
            annual_precip, frtcindx, bgppa, bgppb, agppa, agppb,
            cfrtcw_1, cfrtcw_2, cfrtcn_1, cfrtcn_2)
        self.assert_all_values_in_array_within_range(
            fracrc_p, minimum_acceptable_fracrc_p,
            maximum_acceptable_fracrc_p, _TARGET_NODATA)

        # known values
        annual_precip = numpy.full(array_shape, 42)
        bgppa = numpy.full(array_shape, 101)
        bgppb = numpy.full(array_shape, 4.2)
        agppa = numpy.full(array_shape, -12)
        agppb = numpy.full(array_shape, 3.2)
        cfrtcw_1 = numpy.full(array_shape, 0.4)
        cfrtcw_2 = numpy.full(array_shape, 0.33)
        cfrtcn_1 = numpy.full(array_shape, 0.76)
        cfrtcn_2 = numpy.full(array_shape, 0.02)

        insert_nodata_values_into_array(annual_precip, _TARGET_NODATA)

        known_fracrc_p_frtcindx_0 = 0.69385
        known_fracrc_p_frtcindx_1 = 0.3775
        tolerance = 0.0001

        frtcindx = numpy.full(array_shape, 0)
        fracrc_p = forage.calc_provisional_fracrc(
            annual_precip, frtcindx, bgppa, bgppb, agppa, agppb,
            cfrtcw_1, cfrtcw_2, cfrtcn_1, cfrtcn_2)
        self.assert_all_values_in_array_within_range(
            fracrc_p, known_fracrc_p_frtcindx_0 - tolerance,
            known_fracrc_p_frtcindx_0 + tolerance, _TARGET_NODATA)

        frtcindx = numpy.full(array_shape, 1)
        fracrc_p = forage.calc_provisional_fracrc(
            annual_precip, frtcindx, bgppa, bgppb, agppa, agppb,
            cfrtcw_1, cfrtcw_2, cfrtcn_1, cfrtcn_2)
        self.assert_all_values_in_array_within_range(
            fracrc_p, known_fracrc_p_frtcindx_1 - tolerance,
            known_fracrc_p_frtcindx_1 + tolerance, _TARGET_NODATA)

    def test_calc_ce_ratios(self):
        """Test `calc_ce_ratios`.

        Use the funciton `calc_ce_ratios` to calculate minimum and maximum
        carbon to nutrient ratios. Test that ratios calculated from random
        inputs are within the range of valid values given the range of
        inputs. Introduce nodata values into inputs and test that results
        remain within valid ranges. Calculate the ratios from known inputs
        against results calculated by hand.

        Raises:
            AssertionError if calculated ratios are outside the range of
                valid results given range of random inputs
            AssertionError if calculated ratios are not within 0.0001 of
                values calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        pramn_1_path = os.path.join(self.workspace_dir, 'pramn_1.tif')
        pramn_2_path = os.path.join(self.workspace_dir, 'pramn_2.tif')
        aglivc_path = os.path.join(self.workspace_dir, 'aglivc.tif')
        biomax_path = os.path.join(self.workspace_dir, 'biomax.tif')
        pramx_1_path = os.path.join(self.workspace_dir, 'pramx_1.tif')
        pramx_2_path = os.path.join(self.workspace_dir, 'pramx_2.tif')
        prbmn_1_path = os.path.join(self.workspace_dir, 'prbmn_1.tif')
        prbmn_2_path = os.path.join(self.workspace_dir, 'prbmn_2.tif')
        prbmx_1_path = os.path.join(self.workspace_dir, 'prbmx_1.tif')
        prbmx_2_path = os.path.join(self.workspace_dir, 'prbmx_2.tif')
        annual_precip_path = os.path.join(
            self.workspace_dir, 'annual_precip.tif')
        create_random_raster(pramn_1_path, 20, 50)
        create_random_raster(pramn_2_path, 52, 70)
        create_random_raster(aglivc_path, 20, 400)
        create_random_raster(biomax_path, 300, 500)
        create_random_raster(pramx_1_path, 51, 100)
        create_random_raster(pramx_2_path, 70, 130)
        create_random_raster(prbmn_1_path, 30, 70)
        create_random_raster(prbmn_2_path, 0, 0.2)
        create_random_raster(prbmx_1_path, 40, 70)
        create_random_raster(prbmx_2_path, 0, 0.4)
        create_random_raster(annual_precip_path, 22, 100)

        pft_i = numpy.random.randint(0, 5)
        iel = numpy.random.randint(1, 3)

        month_reg = {
            'cercrp_min_above_{}_{}'.format(iel, pft_i): os.path.join(
                self.workspace_dir,
                'cercrp_min_above_{}_{}.tif'.format(iel, pft_i)),
            'cercrp_max_above_{}_{}'.format(iel, pft_i): os.path.join(
                self.workspace_dir,
                'cercrp_max_above_{}_{}.tif'.format(iel, pft_i)),
            'cercrp_min_below_{}_{}'.format(iel, pft_i): os.path.join(
                self.workspace_dir,
                'cercrp_min_below_{}_{}.tif'.format(iel, pft_i)),
            'cercrp_max_below_{}_{}'.format(iel, pft_i): os.path.join(
                self.workspace_dir,
                'cercrp_max_below_{}_{}.tif'.format(iel, pft_i)),
        }

        acceptable_range_dict = {
            'cercrp_min_above_{}_{}'.format(iel, pft_i): {
                'minimum_acceptable_value': 25.3333,
                'maximum_acceptable_value': 70.,
            },
            'cercrp_max_above_{}_{}'.format(iel, pft_i): {
                'minimum_acceptable_value': 25.,
                'maximum_acceptable_value': 130.,
            },
            'cercrp_min_below_{}_{}'.format(iel, pft_i): {
                'minimum_acceptable_value': 30.,
                'maximum_acceptable_value': 90.,
            },
            'cercrp_max_below_{}_{}'.format(iel, pft_i): {
                'minimum_acceptable_value': 40.,
                'maximum_acceptable_value': 110.,
            },
        }
        forage.calc_ce_ratios(
            pramn_1_path, pramn_2_path, aglivc_path, biomax_path,
            pramx_1_path, pramx_2_path, prbmn_1_path, prbmn_2_path,
            prbmx_1_path, prbmx_2_path, annual_precip_path, month_reg,
            pft_i, iel)
        for path, ranges in acceptable_range_dict.iteritems():
            self.assert_all_values_in_raster_within_range(
                month_reg[path], ranges['minimum_acceptable_value'],
                ranges['maximum_acceptable_value'], _TARGET_NODATA)

        insert_nodata_values_into_raster(aglivc_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(prbmn_1_path, _IC_NODATA)
        insert_nodata_values_into_raster(annual_precip_path, _TARGET_NODATA)
        forage.calc_ce_ratios(
            pramn_1_path, pramn_2_path, aglivc_path, biomax_path,
            pramx_1_path, pramx_2_path, prbmn_1_path, prbmn_2_path,
            prbmx_1_path, prbmx_2_path, annual_precip_path, month_reg,
            pft_i, iel)
        for path, ranges in acceptable_range_dict.iteritems():
            self.assert_all_values_in_raster_within_range(
                month_reg[path], ranges['minimum_acceptable_value'],
                ranges['maximum_acceptable_value'], _TARGET_NODATA)

        # known inputs
        create_random_raster(pramn_1_path, 22, 22)
        create_random_raster(pramn_2_path, 55, 55)
        create_random_raster(aglivc_path, 321, 321)
        create_random_raster(biomax_path, 300, 300)
        create_random_raster(pramx_1_path, 46, 46)
        create_random_raster(pramx_2_path, 78, 78)
        create_random_raster(prbmn_1_path, 52, 52)
        create_random_raster(prbmn_2_path, 0.18, 0.18)
        create_random_raster(prbmx_1_path, 42, 42)
        create_random_raster(prbmx_2_path, 0.33, 0.33)
        create_random_raster(annual_precip_path, 77.22, 77.22)

        known_value_dict = {
            'cercrp_min_above_{}_{}'.format(iel, pft_i): 55.,
            'cercrp_max_above_{}_{}'.format(iel, pft_i): 78.,
            'cercrp_min_below_{}_{}'.format(iel, pft_i): 65.8996,
            'cercrp_max_below_{}_{}'.format(iel, pft_i): 67.4826,
        }
        tolerance = 0.0001

        insert_nodata_values_into_raster(aglivc_path, _SV_NODATA)
        insert_nodata_values_into_raster(prbmn_1_path, _IC_NODATA)
        insert_nodata_values_into_raster(annual_precip_path, _TARGET_NODATA)
        forage.calc_ce_ratios(
            pramn_1_path, pramn_2_path, aglivc_path, biomax_path,
            pramx_1_path, pramx_2_path, prbmn_1_path, prbmn_2_path,
            prbmx_1_path, prbmx_2_path, annual_precip_path, month_reg,
            pft_i, iel)
        for path, value in known_value_dict.iteritems():
            self.assert_all_values_in_raster_within_range(
                month_reg[path], value - tolerance,
                value + tolerance, _TARGET_NODATA)

    def test_calc_revised_fracrc(self):
        """Test `calc_revised_fracrc`.

        Use the function `calc_revised_fracrc` to calculate fracrc_r, the
        revised fraction of carbon allocated to roots. Test that fracrc_r
        calculated from random inputs is within the range of valid values
        according to the range of inputs. Introduce nodata values into
        inputs and test that fracrc_r remains within the valid range.
        Test fracrc_r calculated from known inputs against the result
        calculated by hand.

        Raises:
            AssertionError if fracrc_r calculated from random inputs is
                outside the range of valid values
            AssertionError if fracrc_r from known inputs is not within
                0.0001 of the value calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        frtcindx_path = os.path.join(self.workspace_dir, 'frtcindx.tif')
        fracrc_p_path = os.path.join(self.workspace_dir, 'fracrc_p.tif')
        totale_1_path = os.path.join(self.workspace_dir, 'totale_1.tif')
        totale_2_path = os.path.join(self.workspace_dir, 'totale_2.tif')
        demand_1_path = os.path.join(self.workspace_dir, 'demand_1.tif')
        demand_2_path = os.path.join(self.workspace_dir, 'demand_2.tif')
        h2ogef_1_path = os.path.join(self.workspace_dir, 'h2ogef_1.tif')
        cfrtcw_1_path = os.path.join(self.workspace_dir, 'cfrtcw_1.tif')
        cfrtcw_2_path = os.path.join(self.workspace_dir, 'cfrtcw_2.tif')
        cfrtcn_1_path = os.path.join(self.workspace_dir, 'cfrtcn_1.tif')
        cfrtcn_2_path = os.path.join(self.workspace_dir, 'cfrtcn_2.tif')
        fracrc_r_path = os.path.join(self.workspace_dir, 'fracrc_r.tif')

        create_random_raster(fracrc_p_path, 0.2, 0.95)
        create_random_raster(totale_1_path, 0, 320)
        create_random_raster(totale_2_path, 0, 52)
        create_random_raster(demand_1_path, 1, 34)
        create_random_raster(demand_2_path, 1, 34)
        create_random_raster(h2ogef_1_path, 0.01, 0.9)
        create_random_raster(cfrtcw_1_path, 0.4, 0.8)
        create_random_raster(cfrtcw_2_path, 1.01, 0.39)
        create_random_raster(cfrtcn_1_path, 0.4, 0.8)
        create_random_raster(cfrtcn_2_path, 1.01, 0.39)

        minimum_acceptable_fracrc_r = 0.01
        maximum_acceptable_fracrc_r = 0.999

        create_random_raster(frtcindx_path, 0, 0)
        forage.calc_revised_fracrc(
            frtcindx_path, fracrc_p_path, totale_1_path, totale_2_path,
            demand_1_path, demand_2_path, h2ogef_1_path, cfrtcw_1_path,
            cfrtcw_2_path, cfrtcn_1_path, cfrtcn_2_path, fracrc_r_path)
        self.assert_all_values_in_raster_within_range(
            fracrc_r_path, minimum_acceptable_fracrc_r,
            maximum_acceptable_fracrc_r, _TARGET_NODATA)

        create_random_raster(frtcindx_path, 1, 1)
        forage.calc_revised_fracrc(
            frtcindx_path, fracrc_p_path, totale_1_path, totale_2_path,
            demand_1_path, demand_2_path, h2ogef_1_path, cfrtcw_1_path,
            cfrtcw_2_path, cfrtcn_1_path, cfrtcn_2_path, fracrc_r_path)
        self.assert_all_values_in_raster_within_range(
            fracrc_r_path, minimum_acceptable_fracrc_r,
            maximum_acceptable_fracrc_r, _TARGET_NODATA)

        insert_nodata_values_into_raster(fracrc_p_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(totale_2_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(demand_1_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(cfrtcw_2_path, _IC_NODATA)

        create_random_raster(frtcindx_path, 0, 0)
        forage.calc_revised_fracrc(
            frtcindx_path, fracrc_p_path, totale_1_path, totale_2_path,
            demand_1_path, demand_2_path, h2ogef_1_path, cfrtcw_1_path,
            cfrtcw_2_path, cfrtcn_1_path, cfrtcn_2_path, fracrc_r_path)
        self.assert_all_values_in_raster_within_range(
            fracrc_r_path, minimum_acceptable_fracrc_r,
            maximum_acceptable_fracrc_r, _TARGET_NODATA)

        create_random_raster(frtcindx_path, 1, 1)
        forage.calc_revised_fracrc(
            frtcindx_path, fracrc_p_path, totale_1_path, totale_2_path,
            demand_1_path, demand_2_path, h2ogef_1_path, cfrtcw_1_path,
            cfrtcw_2_path, cfrtcn_1_path, cfrtcn_2_path, fracrc_r_path)
        self.assert_all_values_in_raster_within_range(
            fracrc_r_path, minimum_acceptable_fracrc_r,
            maximum_acceptable_fracrc_r, _TARGET_NODATA)

        # known values
        create_random_raster(fracrc_p_path, 0.7, 0.7)
        create_random_raster(totale_1_path, 10, 10)
        create_random_raster(totale_2_path, 157, 157)
        create_random_raster(demand_1_path, 14, 14)
        create_random_raster(demand_2_path, 30, 30)
        create_random_raster(h2ogef_1_path, 0.7, 0.7)
        create_random_raster(cfrtcw_1_path, 0.7, 0.7)
        create_random_raster(cfrtcw_2_path, 0.36, 0.36)
        create_random_raster(cfrtcn_1_path, 0.47, 0.47)
        create_random_raster(cfrtcn_2_path, 0.33, 0.33)

        known_fracrc_r_frtcindx_0 = 0.7
        known_fracrc_r_frtcindx_1 = 0.462
        tolerance = 0.0001

        create_random_raster(frtcindx_path, 0, 0)
        forage.calc_revised_fracrc(
            frtcindx_path, fracrc_p_path, totale_1_path, totale_2_path,
            demand_1_path, demand_2_path, h2ogef_1_path, cfrtcw_1_path,
            cfrtcw_2_path, cfrtcn_1_path, cfrtcn_2_path, fracrc_r_path)
        self.assert_all_values_in_raster_within_range(
            fracrc_r_path, known_fracrc_r_frtcindx_0 - tolerance,
            known_fracrc_r_frtcindx_0 + tolerance, _TARGET_NODATA)

        create_random_raster(frtcindx_path, 1, 1)
        forage.calc_revised_fracrc(
            frtcindx_path, fracrc_p_path, totale_1_path, totale_2_path,
            demand_1_path, demand_2_path, h2ogef_1_path, cfrtcw_1_path,
            cfrtcw_2_path, cfrtcn_1_path, cfrtcn_2_path, fracrc_r_path)
        self.assert_all_values_in_raster_within_range(
            fracrc_r_path, known_fracrc_r_frtcindx_1 - tolerance,
            known_fracrc_r_frtcindx_1 + tolerance, _TARGET_NODATA)

    def test_grazing_effect(self):
        """Test `grazing_effect_on_aboveground_production`.

        Use the function `grazing_effect_on_aboveground_production` to
        calculate agprod, revised aboveground production including the
        effects of grazing. Use the function `grazing_effect_on_root_shoot` to
        calculate rtsh, root:shoot ratio including the effects of grazing.
        Test that agprod and rtsh match values calculated by hand. Introduce
        nodata values into inputs and ensure that calculated agprod and rtsh
        still match value calculated by hand.

        Raises:
            AssertionError if agprod is not within 0.0001 of value
                calculated by hand
            AssertionError if rtsh if not within 0.0001 of value calculated by
                hand

        Returns:
            None
        """
        from natcap.invest import forage

        array_shape = (3, 3)

        # known values
        tgprod = numpy.full(array_shape, 500)
        fracrc = numpy.full(array_shape, 0.62)
        flgrem = numpy.full(array_shape, 0.16)
        gremb = numpy.full(array_shape, 0.02)

        tolerance = 0.0001

        grzeff = numpy.full(array_shape, 1)
        agprod_grzeff_1 = 122.816
        rtsh_grzeff_1 = 1.63158
        agprod = forage.grazing_effect_on_aboveground_production(
            tgprod, fracrc, flgrem, grzeff)
        rtsh = forage.grazing_effect_on_root_shoot(
            fracrc, flgrem, grzeff, gremb)
        self.assert_all_values_in_array_within_range(
            agprod, agprod_grzeff_1 - tolerance, agprod_grzeff_1 + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            rtsh, rtsh_grzeff_1 - tolerance, rtsh_grzeff_1 + tolerance,
            _TARGET_NODATA)

        grzeff = numpy.full(array_shape, 2)
        agprod_grzeff_2 = 240.6828
        rtsh_grzeff_2 = 1.818
        agprod = forage.grazing_effect_on_aboveground_production(
            tgprod, fracrc, flgrem, grzeff)
        rtsh = forage.grazing_effect_on_root_shoot(
            fracrc, flgrem, grzeff, gremb)
        self.assert_all_values_in_array_within_range(
            agprod, agprod_grzeff_2 - tolerance, agprod_grzeff_2 + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            rtsh, rtsh_grzeff_2 - tolerance, rtsh_grzeff_2 + tolerance,
            _TARGET_NODATA)

        grzeff = numpy.full(array_shape, 3)
        agprod_grzeff_3 = 190
        rtsh_grzeff_3 = 1.818
        agprod = forage.grazing_effect_on_aboveground_production(
            tgprod, fracrc, flgrem, grzeff)
        rtsh = forage.grazing_effect_on_root_shoot(
            fracrc, flgrem, grzeff, gremb)
        self.assert_all_values_in_array_within_range(
            agprod, agprod_grzeff_3 - tolerance, agprod_grzeff_3 + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            rtsh, rtsh_grzeff_3 - tolerance, rtsh_grzeff_3 + tolerance,
            _TARGET_NODATA)

        grzeff = numpy.full(array_shape, 4)
        agprod_grzeff_4 = 190
        rtsh_grzeff_4 = 0.9968
        agprod = forage.grazing_effect_on_aboveground_production(
            tgprod, fracrc, flgrem, grzeff)
        rtsh = forage.grazing_effect_on_root_shoot(
            fracrc, flgrem, grzeff, gremb)
        self.assert_all_values_in_array_within_range(
            agprod, agprod_grzeff_4 - tolerance, agprod_grzeff_4 + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            rtsh, rtsh_grzeff_4 - tolerance, rtsh_grzeff_4 + tolerance,
            _TARGET_NODATA)

        grzeff = numpy.full(array_shape, 5)
        agprod_grzeff_5 = 240.6828
        rtsh_grzeff_5 = 0.9968
        agprod = forage.grazing_effect_on_aboveground_production(
            tgprod, fracrc, flgrem, grzeff)
        rtsh = forage.grazing_effect_on_root_shoot(
            fracrc, flgrem, grzeff, gremb)
        self.assert_all_values_in_array_within_range(
            agprod, agprod_grzeff_5 - tolerance, agprod_grzeff_5 + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            rtsh, rtsh_grzeff_5 - tolerance, rtsh_grzeff_5 + tolerance,
            _TARGET_NODATA)

        grzeff = numpy.full(array_shape, 6)
        agprod_grzeff_6 = 122.816
        rtsh_grzeff_6 = 0.9968
        agprod = forage.grazing_effect_on_aboveground_production(
            tgprod, fracrc, flgrem, grzeff)
        rtsh = forage.grazing_effect_on_root_shoot(
            fracrc, flgrem, grzeff, gremb)
        self.assert_all_values_in_array_within_range(
            agprod, agprod_grzeff_6 - tolerance, agprod_grzeff_6 + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            rtsh, rtsh_grzeff_6 - tolerance, rtsh_grzeff_6 + tolerance,
            _TARGET_NODATA)

        insert_nodata_values_into_array(fracrc, _TARGET_NODATA)

        grzeff = numpy.full(array_shape, 4)
        agprod_grzeff_4 = 190
        rtsh_grzeff_4 = 0.9968
        agprod = forage.grazing_effect_on_aboveground_production(
            tgprod, fracrc, flgrem, grzeff)
        rtsh = forage.grazing_effect_on_root_shoot(
            fracrc, flgrem, grzeff, gremb)
        self.assert_all_values_in_array_within_range(
            agprod, agprod_grzeff_4 - tolerance, agprod_grzeff_4 + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            rtsh, rtsh_grzeff_4 - tolerance, rtsh_grzeff_4 + tolerance,
            _TARGET_NODATA)

        grzeff = numpy.full(array_shape, 2)
        agprod_grzeff_2 = 240.6828
        rtsh_grzeff_2 = 1.818
        agprod = forage.grazing_effect_on_aboveground_production(
            tgprod, fracrc, flgrem, grzeff)
        rtsh = forage.grazing_effect_on_root_shoot(
            fracrc, flgrem, grzeff, gremb)
        self.assert_all_values_in_array_within_range(
            agprod, agprod_grzeff_2 - tolerance, agprod_grzeff_2 + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            rtsh, rtsh_grzeff_2 - tolerance, rtsh_grzeff_2 + tolerance,
            _TARGET_NODATA)

    def test_calc_tgprod_final(self):
        """Test `calc_tgprod_final`.

        Use the function `calc_tgprod_final` to calculate tgprod, final total
        prodcution from root:shoot ratio and aboveground production. Test that
        calculated tgprod matches results calculated by hand.

        Raises:
            AssertionError if tgprod is not within 0.0001 of the value
                calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        array_size = (3, 3)
        # known values
        rtsh = numpy.full(array_size, 0.72)
        agprod = numpy.full(array_size, 333)

        known_tgprod = 572.76
        tolerance = 0.0001
        tgprod = forage.calc_tgprod_final(rtsh, agprod)
        self.assert_all_values_in_array_within_range(
            tgprod, known_tgprod - tolerance, known_tgprod + tolerance,
            _TARGET_NODATA)

    def test_snow(self):
        """Test `_snow`.

        Use the function `_snow` to modify snow pack, evaporate from snow pack,
        melt snow, and determine liquid inputs to soil after snow. Test
        the raster-based function against a point-based function defined
        here.

        Raises:
            AssertionError if raster-based outputs do not match outputs
                calculated by point-based version

        Returns:
            None
        """
        def snow_point(
                precip, max_temp, min_temp, snow, snlq, pet, tmelt_1, tmelt_2,
                shwave):
            """Point-based implementation of `_snow`.

            This implementation reproduces Century's process for determining
            snowfall, evaporation from snow, snowmelt, and liquid draining
            into soil after snow is accounted.

            Parameters:
                precip (float): precipitation this month
                max_temp (float): maximum temperature this month
                min_temp (float): minimum temperature this month
                snow (float): existing snow prior to new precipitation
                snlq (float): existing liquid in snow prior to new
                    precipitation
                pet (float): potential evapotranspiration
                tmelt_1 (float): parameter, temperature above which some
                    snow melts
                tmelt_2 (float): parameter, ratio between degrees above the
                    minimum temperature and cm of snow that will melt

            Returns:
                dict of modified quantities: snowmelt, snow, snlq, pet,
                inputs_after_snow
            """
            tave = (max_temp + min_temp) / 2.
            inputs_after_snow = precip
            # precip falls as snow when temperature is below freezing
            if tave <= 0:
                snow = snow + precip
                # all precip is snow, none left
                inputs_after_snow = 0
            else:
                if snow > 0:
                    snlq = snlq + precip
                    # all precip is rain on snow, none left
                    inputs_after_snow = 0

            snowmelt = 0
            if snow > 0:
                snowtot = snow + snlq
                evsnow = min(snowtot, (pet * 0.87))
                snow = max(snow - evsnow * (snow/snowtot), 0.)
                snlq = max(snlq-evsnow * (snlq/snowtot), 0.)
                pet = max((pet - evsnow / 0.87), 0)

                if tave >= tmelt_1:
                    snowmelt = tmelt_2 * (tave - tmelt_1) * shwave
                    snowmelt = max(snowmelt, 0)
                    snowmelt = min(snowmelt, snow)
                    snow = snow - snowmelt
                    snlq = snlq + snowmelt
                    if snlq > (0.5 * snow):
                        add = snlq - 0.5 * snow
                        snlq = snlq - add
                        inputs_after_snow = add

            results_dict = {
                'snowmelt': snowmelt,
                'snow': snow,
                'snlq': snlq,
                'pet': pet,
                'inputs_after_snow': inputs_after_snow,
            }
            return results_dict

        from natcap.invest import forage

        # shortwave radiation and pet calculated by hand
        CURRENT_MONTH = 10
        SHWAVE = 437.04
        TMELT_1 = 0.
        TMELT_2 = 0.002
        FWLOSS_4 = 0.6

        # rain on snow, all snow melts
        test_dict = {
            'precip': 15.,
            'max_temp': 23.,
            'min_temp': -2.,
            'snow': 8.,
            'snlq': 4.,
            'pet': 4.967985,
            'tmelt_1': TMELT_1,
            'tmelt_2': TMELT_2,
            'shwave': SHWAVE,
        }
        test_dict['tave'] = (test_dict['max_temp'] + test_dict['min_temp']) / 2

        site_param_table = {
            1: {
                'tmelt_1': TMELT_1,
                'tmelt_2': TMELT_2,
                'fwloss_4': FWLOSS_4,
            }
        }
        site_index_path = os.path.join(self.workspace_dir, 'site_index.tif')
        precip_path = os.path.join(self.workspace_dir, 'precip.tif')
        tave_path = os.path.join(self.workspace_dir, 'tave.tif')
        max_temp_path = os.path.join(self.workspace_dir, 'max_temp.tif')
        min_temp_path = os.path.join(self.workspace_dir, 'min_temp.tif')
        prev_snow_path = os.path.join(self.workspace_dir, 'prev_snow.tif')
        prev_snlq_path = os.path.join(self.workspace_dir, 'prev_snlq.tif')
        snowmelt_path = os.path.join(self.workspace_dir, 'snowmelt.tif')
        snow_path = os.path.join(self.workspace_dir, 'snow.tif')
        snlq_path = os.path.join(self.workspace_dir, 'snlq.tif')
        inputs_after_snow_path = os.path.join(
            self.workspace_dir, 'inputs_after_snow.tif')
        pet_rem_path = os.path.join(self.workspace_dir, 'pet_rem.tif')

        # raster inputs
        nrows = 1
        ncols = 1
        create_random_raster(site_index_path, 1, 1, nrows=nrows, ncols=ncols)
        create_random_raster(
            precip_path, test_dict['precip'], test_dict['precip'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            tave_path, test_dict['tave'], test_dict['tave'], nrows=nrows,
            ncols=ncols)
        create_random_raster(
            max_temp_path, test_dict['max_temp'], test_dict['max_temp'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            min_temp_path, test_dict['min_temp'], test_dict['min_temp'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            prev_snow_path, test_dict['snow'], test_dict['snow'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            prev_snlq_path, test_dict['snlq'], test_dict['snlq'],
            nrows=nrows, ncols=ncols)

        tolerance = 0.000015
        result_dict = snow_point(
            test_dict['precip'], test_dict['max_temp'], test_dict['min_temp'],
            test_dict['snow'], test_dict['snlq'], test_dict['pet'],
            test_dict['tmelt_1'], test_dict['tmelt_2'], SHWAVE)

        forage._snow(
            site_index_path, site_param_table, precip_path, tave_path,
            max_temp_path, min_temp_path, prev_snow_path, prev_snlq_path,
            CURRENT_MONTH, snowmelt_path, snow_path, snlq_path,
            inputs_after_snow_path, pet_rem_path)

        self.assert_all_values_in_raster_within_range(
            pet_rem_path, result_dict['pet'] - tolerance,
            result_dict['pet'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snowmelt_path, result_dict['snowmelt'] - tolerance,
            result_dict['snowmelt'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snow_path, result_dict['snow'], result_dict['snow'],
            _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snlq_path, result_dict['snlq'], result_dict['snlq'],
            _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            inputs_after_snow_path,
            result_dict['inputs_after_snow'] - tolerance,
            result_dict['inputs_after_snow'] + tolerance, _TARGET_NODATA)

        # new snowfall, no snowmelt
        test_dict = {
            'precip': 15.,
            'max_temp': 2.,
            'min_temp': -6.,
            'snow': 2.,
            'snlq': 1.,
            'pet': 1.5690107,
            'tmelt_1': TMELT_1,
            'tmelt_2': TMELT_2,
            'shwave': SHWAVE,
        }
        test_dict['tave'] = (test_dict['max_temp'] + test_dict['min_temp']) / 2

        create_random_raster(site_index_path, 1, 1, nrows=nrows, ncols=ncols)
        create_random_raster(
            precip_path, test_dict['precip'], test_dict['precip'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            tave_path, test_dict['tave'], test_dict['tave'], nrows=nrows,
            ncols=ncols)
        create_random_raster(
            max_temp_path, test_dict['max_temp'], test_dict['max_temp'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            min_temp_path, test_dict['min_temp'], test_dict['min_temp'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            prev_snow_path, test_dict['snow'], test_dict['snow'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            prev_snlq_path, test_dict['snlq'], test_dict['snlq'],
            nrows=nrows, ncols=ncols)

        result_dict = snow_point(
            test_dict['precip'], test_dict['max_temp'], test_dict['min_temp'],
            test_dict['snow'], test_dict['snlq'], test_dict['pet'],
            test_dict['tmelt_1'], test_dict['tmelt_2'], SHWAVE)

        forage._snow(
            site_index_path, site_param_table, precip_path, tave_path,
            max_temp_path, min_temp_path, prev_snow_path, prev_snlq_path,
            CURRENT_MONTH, snowmelt_path, snow_path, snlq_path,
            inputs_after_snow_path, pet_rem_path)

        self.assert_all_values_in_raster_within_range(
            snowmelt_path, result_dict['snowmelt'] - tolerance,
            result_dict['snowmelt'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snow_path, result_dict['snow'] - tolerance,
            result_dict['snow'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snlq_path, result_dict['snlq'] - tolerance,
            result_dict['snlq'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            inputs_after_snow_path, result_dict['inputs_after_snow'],
            result_dict['inputs_after_snow'], _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            pet_rem_path, result_dict['pet'], result_dict['pet'],
            _TARGET_NODATA)

        # new snowfall, some snowmelt
        test_dict = {
            'precip': 10.,
            'max_temp': 2.,
            'min_temp': -2.01,
            'snow': 2.,
            'snlq': 1.,
            'pet': 1.2511096,
            'tmelt_1': TMELT_1,
            'tmelt_2': TMELT_2,
            'shwave': SHWAVE,
        }
        test_dict['tave'] = (test_dict['max_temp'] + test_dict['min_temp']) / 2

        create_random_raster(site_index_path, 1, 1, nrows=nrows, ncols=ncols)
        create_random_raster(
            precip_path, test_dict['precip'], test_dict['precip'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            tave_path, test_dict['tave'], test_dict['tave'], nrows=nrows,
            ncols=ncols)
        create_random_raster(
            max_temp_path, test_dict['max_temp'], test_dict['max_temp'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            min_temp_path, test_dict['min_temp'], test_dict['min_temp'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            prev_snow_path, test_dict['snow'], test_dict['snow'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            prev_snlq_path, test_dict['snlq'], test_dict['snlq'],
            nrows=nrows, ncols=ncols)

        result_dict = snow_point(
            test_dict['precip'], test_dict['max_temp'], test_dict['min_temp'],
            test_dict['snow'], test_dict['snlq'], test_dict['pet'],
            test_dict['tmelt_1'], test_dict['tmelt_2'], SHWAVE)

        forage._snow(
            site_index_path, site_param_table, precip_path, tave_path,
            max_temp_path, min_temp_path, prev_snow_path, prev_snlq_path,
            CURRENT_MONTH, snowmelt_path, snow_path, snlq_path,
            inputs_after_snow_path, pet_rem_path)

        self.assert_all_values_in_raster_within_range(
            snowmelt_path, result_dict['snowmelt'] - tolerance,
            result_dict['snowmelt'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snow_path, result_dict['snow'] - tolerance,
            result_dict['snow'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snlq_path, result_dict['snlq'] - tolerance,
            result_dict['snlq'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            inputs_after_snow_path,
            result_dict['inputs_after_snow'],
            result_dict['inputs_after_snow'], _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            pet_rem_path, result_dict['pet'] - tolerance,
            result_dict['pet'] + tolerance, _TARGET_NODATA)

        # no snow, no snowfall
        test_dict = {
            'precip': 10.,
            'max_temp': 23.,
            'min_temp': -2.,
            'snow': 0.,
            'snlq': 0.,
            'pet': 4.9680004,
            'tmelt_1': TMELT_1,
            'tmelt_2': TMELT_2,
            'shwave': 457.95056,
        }
        test_dict['tave'] = (test_dict['max_temp'] + test_dict['min_temp']) / 2

        create_random_raster(site_index_path, 1, 1, nrows=nrows, ncols=ncols)
        create_random_raster(
            precip_path, test_dict['precip'], test_dict['precip'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            tave_path, test_dict['tave'], test_dict['tave'], nrows=nrows,
            ncols=ncols)
        create_random_raster(
            max_temp_path, test_dict['max_temp'], test_dict['max_temp'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            min_temp_path, test_dict['min_temp'], test_dict['min_temp'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            prev_snow_path, test_dict['snow'], test_dict['snow'],
            nrows=nrows, ncols=ncols)
        create_random_raster(
            prev_snlq_path, test_dict['snlq'], test_dict['snlq'],
            nrows=nrows, ncols=ncols)

        result_dict = snow_point(
            test_dict['precip'], test_dict['max_temp'], test_dict['min_temp'],
            test_dict['snow'], test_dict['snlq'], test_dict['pet'],
            test_dict['tmelt_1'], test_dict['tmelt_2'], SHWAVE)

        forage._snow(
            site_index_path, site_param_table, precip_path, tave_path,
            max_temp_path, min_temp_path, prev_snow_path, prev_snlq_path,
            CURRENT_MONTH, snowmelt_path, snow_path, snlq_path,
            inputs_after_snow_path, pet_rem_path)

        self.assert_all_values_in_raster_within_range(
            snowmelt_path, result_dict['snowmelt'] - tolerance,
            result_dict['snowmelt'] + tolerance, _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snow_path, result_dict['snow'], result_dict['snow'],
            _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            snlq_path, result_dict['snlq'], result_dict['snlq'],
            _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            inputs_after_snow_path, result_dict['inputs_after_snow'],
            result_dict['inputs_after_snow'], _TARGET_NODATA)

        self.assert_all_values_in_raster_within_range(
            pet_rem_path, result_dict['pet'] - tolerance,
            result_dict['pet'] + tolerance, _TARGET_NODATA)

    def test_calc_aboveground_live_biomass(self):
        """Test `_calc_aboveground_live_biomass`.

        Use the function `_calc_aboveground_live_biomass` to calculate
        aboveground live biomass for the purposes of soil water. Test that the
        function reproduces results calculated by hand.

        Raises:
            AssertionError if `_calc_aboveground_live_biomass` does not match
                results calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        array_size = (3, 3)
        # known values
        sum_aglivc = numpy.full(array_size, 200.)
        sum_tgprod = numpy.full(array_size, 180.)

        known_aliv = 545.
        tolerance = 0.00001
        aliv = forage._calc_aboveground_live_biomass(sum_aglivc, sum_tgprod)
        self.assert_all_values_in_array_within_range(
            aliv, known_aliv - tolerance, known_aliv + tolerance,
            _TARGET_NODATA)

    def test_calc_standing_biomass(self):
        """Test `_calc_standing_biomass`.

        Use the function `_calc_standing_biomass` to calculate total
        aboveground standing biomass for soil water. Test that the function
        reproduces results calculated by hand.

        Raises:
            AssertionError if `_calc_standing_biomass` does not match results
                calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        array_size = (3, 3)
        # known values
        aliv = numpy.full(array_size, 545)
        sum_stdedc = numpy.full(array_size, 232)

        known_sd = 800.
        tolerance = 0.00001
        sd = forage._calc_standing_biomass(aliv, sum_stdedc)
        self.assert_all_values_in_array_within_range(
            sd, known_sd - tolerance, known_sd + tolerance, _TARGET_NODATA)

        # known values
        aliv = numpy.full(array_size, 233.2)
        sum_stdedc = numpy.full(array_size, 172)

        known_sd = 663.2
        tolerance = 0.0001
        sd = forage._calc_standing_biomass(aliv, sum_stdedc)
        self.assert_all_values_in_array_within_range(
            sd, known_sd - tolerance, known_sd + tolerance, _TARGET_NODATA)

    def test_subtract_surface_losses(self):
        """Test `subtract_surface_losses`.

        Use the function `subtract_surface_losses` to calculate moisture
        losses to runoff, canopy interception, and evaporation.  Test that
        the function reproduces results calculated by a point-based function
        defined here.

        Raises:
            AssertionError if `subtract_surface_losses` does not match
                point-based results calculated by test function

        Returns:
            None
        """
        def surface_losses_point(
                inputs_after_snow, fracro, precro, snow, alit, sd, fwloss_1,
                fwloss_2, pet_rem):
            """Point- based implementation of `subtract_surface_losses`.

            This implementation reproduces Century's process for determining
            loss of moisture inputs to runoff, canopy interception, and surface
            evaporation.

            Parameters:
                inputs_after_snow (float): surface water inputs from
                    precipitation and snowmelt, prior to runoff
                fracro (float): parameter, fraction of surface water
                    above precro that is lost to runoff
                precro (float): parameter, amount of surface water that
                    must be available for runoff to occur
                snow (float): current snowpack
                alit (float): biomass in surface litter
                sd (float): total standing biomass
                fwloss_1 (float): parameter, scaling factor for
                    interception and evaporation of precip by vegetation
                fwloss_2 (float): parameter, scaling factor for bare soil
                    evaporation of precip
                pet_rem (float): potential evaporation remaining after
                    evaporation of snow

            Returns:
                dict of modified quantities: inputs_after_surface, surface
                    water inputs to soil after runoff and surface evaporation
                    are subtracted; absevap, moisture lost to surface
                    evaporation; and evap_losses, total surface evaporation
            """
            runoff = max(fracro * (inputs_after_snow - precro), 0.)
            inputs_after_runoff = inputs_after_snow - runoff
            if snow == 0:
                aint = (0.0003 * alit + 0.0006 * sd) * fwloss_1
                absevap = (
                    0.5 * math.exp((-0.002 * alit) - (0.004 * sd)) * fwloss_2)
                evap_losses = min(
                    ((absevap + aint) * inputs_after_runoff), 0.4 * pet_rem)
            else:
                absevap = 0
                evap_losses = 0
            inputs_after_surface = inputs_after_runoff - evap_losses

            results_dict = {
                'inputs_after_surface': inputs_after_surface,
                'absevap': absevap,
                'evap_losses': evap_losses,
            }
            return results_dict

        from natcap.invest import forage
        array_size = (10, 10)

        # snow cover, runoff losses only
        test_dict = {
            'inputs_after_snow': 34.,
            'fracro': 0.15,
            'precro': 8.,
            'snow': 20.,
            'alit': 100.,
            'sd': 202.5,
            'fwloss_1': 0.9,
            'fwloss_2': 0.7,
            'pet_rem': 3.88,
        }

        inputs_after_snow = numpy.full(
            array_size, test_dict['inputs_after_snow'])
        fracro = numpy.full(array_size, test_dict['fracro'])
        precro = numpy.full(array_size, test_dict['precro'])
        snow = numpy.full(array_size, test_dict['snow'])
        alit = numpy.full(array_size, test_dict['alit'])
        sd = numpy.full(array_size, test_dict['sd'])
        fwloss_1 = numpy.full(array_size, test_dict['fwloss_1'])
        fwloss_2 = numpy.full(array_size, test_dict['fwloss_2'])
        pet_rem = numpy.full(array_size, test_dict['pet_rem'])

        result_dict = surface_losses_point(
            test_dict['inputs_after_snow'], test_dict['fracro'],
            test_dict['precro'], test_dict['snow'], test_dict['alit'],
            test_dict['sd'], test_dict['fwloss_1'], test_dict['fwloss_2'],
            test_dict['pet_rem'])
        tolerance = 0.00001
        inputs_after_surface = forage.subtract_surface_losses(
            'inputs_after_surface')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        absevap = forage.subtract_surface_losses(
            'absevap')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        evap_losses = forage.subtract_surface_losses(
            'evap_losses')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)

        self.assert_all_values_in_array_within_range(
            inputs_after_surface,
            result_dict['inputs_after_surface'] - tolerance,
            result_dict['inputs_after_surface'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            absevap, result_dict['absevap'] - tolerance,
            result_dict['absevap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            evap_losses, result_dict['evap_losses'] - tolerance,
            result_dict['evap_losses'] + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(pet_rem, _TARGET_NODATA)
        insert_nodata_values_into_array(snow, _TARGET_NODATA)
        insert_nodata_values_into_array(fracro, _IC_NODATA)

        inputs_after_surface = forage.subtract_surface_losses(
            'inputs_after_surface')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        absevap = forage.subtract_surface_losses(
            'absevap')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        evap_losses = forage.subtract_surface_losses(
            'evap_losses')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)

        self.assert_all_values_in_array_within_range(
            inputs_after_surface,
            result_dict['inputs_after_surface'] - tolerance,
            result_dict['inputs_after_surface'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            absevap, result_dict['absevap'] - tolerance,
            result_dict['absevap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            evap_losses, result_dict['evap_losses'] - tolerance,
            result_dict['evap_losses'] + tolerance, _TARGET_NODATA)

        # no snow cover, large surface biomass
        test_dict = {
            'inputs_after_snow': 12.,
            'fracro': 0.15,
            'precro': 8.,
            'snow': 0.,
            'alit': 200.1,
            'sd': 800.,
            'fwloss_1': 0.8,
            'fwloss_2': 0.8,
            'pet_rem': 3.88,
        }

        inputs_after_snow = numpy.full(
            array_size, test_dict['inputs_after_snow'])
        fracro = numpy.full(array_size, test_dict['fracro'])
        precro = numpy.full(array_size, test_dict['precro'])
        snow = numpy.full(array_size, test_dict['snow'])
        alit = numpy.full(array_size, test_dict['alit'])
        sd = numpy.full(array_size, test_dict['sd'])
        fwloss_1 = numpy.full(array_size, test_dict['fwloss_1'])
        fwloss_2 = numpy.full(array_size, test_dict['fwloss_2'])
        pet_rem = numpy.full(array_size, test_dict['pet_rem'])

        result_dict = surface_losses_point(
            test_dict['inputs_after_snow'], test_dict['fracro'],
            test_dict['precro'], test_dict['snow'], test_dict['alit'],
            test_dict['sd'], test_dict['fwloss_1'], test_dict['fwloss_2'],
            test_dict['pet_rem'])
        tolerance = 0.00001
        inputs_after_surface = forage.subtract_surface_losses(
            'inputs_after_surface')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        absevap = forage.subtract_surface_losses(
            'absevap')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        evap_losses = forage.subtract_surface_losses(
            'evap_losses')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)

        self.assert_all_values_in_array_within_range(
            inputs_after_surface,
            result_dict['inputs_after_surface'] - tolerance,
            result_dict['inputs_after_surface'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            absevap, result_dict['absevap'] - tolerance,
            result_dict['absevap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            evap_losses, result_dict['evap_losses'] - tolerance,
            result_dict['evap_losses'] + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(alit, _TARGET_NODATA)
        insert_nodata_values_into_array(precro, _IC_NODATA)
        insert_nodata_values_into_array(fwloss_2, _IC_NODATA)

        inputs_after_surface = forage.subtract_surface_losses(
            'inputs_after_surface')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        absevap = forage.subtract_surface_losses(
            'absevap')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        evap_losses = forage.subtract_surface_losses(
            'evap_losses')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)

        self.assert_all_values_in_array_within_range(
            inputs_after_surface,
            result_dict['inputs_after_surface'] - tolerance,
            result_dict['inputs_after_surface'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            absevap, result_dict['absevap'] - tolerance,
            result_dict['absevap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            evap_losses, result_dict['evap_losses'] - tolerance,
            result_dict['evap_losses'] + tolerance, _TARGET_NODATA)

        # no snow cover, small surface biomass
        test_dict = {
            'inputs_after_snow': 12.,
            'fracro': 0.15,
            'precro': 8.,
            'snow': 0.,
            'alit': 300.1,
            'sd': 80.5,
            'fwloss_1': 0.8,
            'fwloss_2': 0.8,
            'pet_rem': 4.99,
        }

        inputs_after_snow = numpy.full(
            array_size, test_dict['inputs_after_snow'])
        fracro = numpy.full(array_size, test_dict['fracro'])
        precro = numpy.full(array_size, test_dict['precro'])
        snow = numpy.full(array_size, test_dict['snow'])
        alit = numpy.full(array_size, test_dict['alit'])
        sd = numpy.full(array_size, test_dict['sd'])
        fwloss_1 = numpy.full(array_size, test_dict['fwloss_1'])
        fwloss_2 = numpy.full(array_size, test_dict['fwloss_2'])
        pet_rem = numpy.full(array_size, test_dict['pet_rem'])

        result_dict = surface_losses_point(
            test_dict['inputs_after_snow'], test_dict['fracro'],
            test_dict['precro'], test_dict['snow'], test_dict['alit'],
            test_dict['sd'], test_dict['fwloss_1'], test_dict['fwloss_2'],
            test_dict['pet_rem'])
        tolerance = 0.00001
        inputs_after_surface = forage.subtract_surface_losses(
            'inputs_after_surface')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        absevap = forage.subtract_surface_losses(
            'absevap')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        evap_losses = forage.subtract_surface_losses(
            'evap_losses')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)

        self.assert_all_values_in_array_within_range(
            inputs_after_surface,
            result_dict['inputs_after_surface'] - tolerance,
            result_dict['inputs_after_surface'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            absevap, result_dict['absevap'] - tolerance,
            result_dict['absevap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            evap_losses, result_dict['evap_losses'] - tolerance,
            result_dict['evap_losses'] + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(inputs_after_snow, _TARGET_NODATA)
        insert_nodata_values_into_array(fwloss_1, _IC_NODATA)
        insert_nodata_values_into_array(fwloss_2, _IC_NODATA)

        inputs_after_surface = forage.subtract_surface_losses(
            'inputs_after_surface')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        absevap = forage.subtract_surface_losses(
            'absevap')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)
        evap_losses = forage.subtract_surface_losses(
            'evap_losses')(
                inputs_after_snow, fracro, precro, snow,
                alit, sd, fwloss_1, fwloss_2, pet_rem)

        self.assert_all_values_in_array_within_range(
            inputs_after_surface,
            result_dict['inputs_after_surface'] - tolerance,
            result_dict['inputs_after_surface'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            absevap, result_dict['absevap'] - tolerance,
            result_dict['absevap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            evap_losses, result_dict['evap_losses'] - tolerance,
            result_dict['evap_losses'] + tolerance, _TARGET_NODATA)

    def test_calc_potential_transpiration(self):
        """Test `calc_potential_transpiration`.

        Use the function `calc_potential_transpiration` to calculate
        potential transpiration from the soil by plants, potential  evaporation
        of soil moisture from soil layer 1, initial transpiration water loss,
        and modified water inputs. Test the function against a point-based
        version defined here.

        Raises:
            AssertionError if point-based test version of
                `subtract_surface_losses` does not match values calculated by
                hand
            AssertionError if `subtract_surface_losses` does not match results
                calculated by point-based test version

        Returns:
            None
        """
        def potential_transpiration_point(
                pet_rem, evap_losses, tave, aliv, current_moisture_inputs):
            """Calculate trap, pevp, and modified moisture inputs.

            Parameters:
                pet_rem (float): potential evapotranspiration remaining after
                    evaporation of snow
                evap_losses (float): total surface evaporation
                tave (float): average temperature
                aliv (float): aboveground live biomass
                current_moisture_inputs (float): moisture inputs after surface
                    losses

            Returns:
                dict of moidified quantities: trap, potential transpiration;
                    pevp, potential evaporation from surface soil layer;
                    modified_moisture_inputs, water to be added to soil layers
                    before transpiration losses are accounted
            """
            trap = pet_rem - evap_losses
            if tave < 2:
                pttr = 0
            else:
                pttr = pet_rem * 0.65 * (1 - math.exp(-0.02 * aliv))
            if pttr <= trap:
                trap = pttr
            if trap <= 0:
                trap = 0.01
            pevp = max(pet_rem - trap - evap_losses, 0.)
            tran = min(trap - 0.01, current_moisture_inputs)
            trap = trap - tran
            modified_moisture_inputs = current_moisture_inputs - tran

            results_dict = {
                'trap': trap,
                'pevp': pevp,
                'modified_moisture_inputs': modified_moisture_inputs,
            }
            return results_dict

        from natcap.invest import forage
        array_size = (10, 10)

        # high transpiration limited by water inputs
        test_dict = {
            'pet_rem': 13.2,
            'evap_losses': 5.28,
            'tave': 22.3,
            'aliv': 100.,
            'current_moisture_inputs': 7.4,
        }

        pet_rem = numpy.full(array_size, test_dict['pet_rem'])
        evap_losses = numpy.full(array_size, test_dict['evap_losses'])
        tave = numpy.full(array_size, test_dict['tave'])
        aliv = numpy.full(array_size, test_dict['aliv'])
        current_moisture_inputs = numpy.full(
            array_size, test_dict['current_moisture_inputs'])

        result_dict = potential_transpiration_point(
            test_dict['pet_rem'], test_dict['evap_losses'], test_dict['tave'],
            test_dict['aliv'], test_dict['current_moisture_inputs'])

        insert_nodata_values_into_array(aliv, _TARGET_NODATA)
        insert_nodata_values_into_array(
            current_moisture_inputs, _TARGET_NODATA)
        insert_nodata_values_into_array(tave, _TARGET_NODATA)

        trap = forage.calc_potential_transpiration(
            'trap')(pet_rem, evap_losses, tave, aliv, current_moisture_inputs)
        pevp = forage.calc_potential_transpiration(
            'pevp')(pet_rem, evap_losses, tave, aliv, current_moisture_inputs)
        modified_moisture_inputs = forage.calc_potential_transpiration(
            'modified_moisture_inputs')(
            pet_rem, evap_losses, tave, aliv, current_moisture_inputs)

        # known values calculated by hand
        known_trap = 0.018823
        known_pevp = 0.50118
        known_modified_moisture_inputs = 0
        tolerance = 0.00001

        self.assertAlmostEqual(
            result_dict['trap'], known_trap, delta=tolerance,
            msg=(
                "trap calculated by point-based test version does not match" +
                " value calculated by hand"))
        self.assertAlmostEqual(
            result_dict['pevp'], known_pevp, delta=tolerance,
            msg=(
                "pevp calculated by point-based test version does not match" +
                " value calculated by hand"))
        self.assertAlmostEqual(
            result_dict['modified_moisture_inputs'],
            known_modified_moisture_inputs,
            delta=tolerance,
            msg=(
                "modified_moisture_inputs calculated by point-based test " +
                "version does not match value calculated by hand"))
        self.assert_all_values_in_array_within_range(
            trap, result_dict['trap'] - tolerance,
            result_dict['trap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            pevp, result_dict['pevp'] - tolerance,
            result_dict['pevp'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            modified_moisture_inputs,
            result_dict['modified_moisture_inputs'] - tolerance,
            result_dict['modified_moisture_inputs'] + tolerance,
            _TARGET_NODATA)

        # low temperature, no transpiration occurs
        test_dict = {
            'pet_rem': 3.24,
            'evap_losses': 1.12,
            'tave': 1.,
            'aliv': 180.,
            'current_moisture_inputs': 62.,
        }

        pet_rem = numpy.full(array_size, test_dict['pet_rem'])
        evap_losses = numpy.full(array_size, test_dict['evap_losses'])
        tave = numpy.full(array_size, test_dict['tave'])
        aliv = numpy.full(array_size, test_dict['aliv'])
        current_moisture_inputs = numpy.full(
            array_size, test_dict['current_moisture_inputs'])

        result_dict = potential_transpiration_point(
            test_dict['pet_rem'], test_dict['evap_losses'], test_dict['tave'],
            test_dict['aliv'], test_dict['current_moisture_inputs'])

        trap = forage.calc_potential_transpiration(
            'trap')(pet_rem, evap_losses, tave, aliv, current_moisture_inputs)
        pevp = forage.calc_potential_transpiration(
            'pevp')(pet_rem, evap_losses, tave, aliv, current_moisture_inputs)
        modified_moisture_inputs = forage.calc_potential_transpiration(
            'modified_moisture_inputs')(
            pet_rem, evap_losses, tave, aliv, current_moisture_inputs)

        self.assert_all_values_in_array_within_range(
            trap, result_dict['trap'] - tolerance,
            result_dict['trap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            pevp, result_dict['pevp'] - tolerance,
            result_dict['pevp'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            modified_moisture_inputs,
            result_dict['modified_moisture_inputs'] - tolerance,
            result_dict['modified_moisture_inputs'] + tolerance,
            _TARGET_NODATA)

        insert_nodata_values_into_array(pet_rem, _TARGET_NODATA)
        insert_nodata_values_into_array(evap_losses, _TARGET_NODATA)
        insert_nodata_values_into_array(tave, _TARGET_NODATA)

        trap = forage.calc_potential_transpiration(
            'trap')(pet_rem, evap_losses, tave, aliv, current_moisture_inputs)
        pevp = forage.calc_potential_transpiration(
            'pevp')(pet_rem, evap_losses, tave, aliv, current_moisture_inputs)
        modified_moisture_inputs = forage.calc_potential_transpiration(
            'modified_moisture_inputs')(
            pet_rem, evap_losses, tave, aliv, current_moisture_inputs)

        self.assert_all_values_in_array_within_range(
            trap, result_dict['trap'] - tolerance,
            result_dict['trap'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            pevp, result_dict['pevp'] - tolerance,
            result_dict['pevp'] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            modified_moisture_inputs,
            result_dict['modified_moisture_inputs'] - tolerance,
            result_dict['modified_moisture_inputs'] + tolerance,
            _TARGET_NODATA)

    def test_get_nlaypg_max(self):
        """Test `_get_nlaypg_max`.

        Use the function `get_nlaypg_max` to retrieve the maximum value of
        nlaypg across plant functional types.

        Raises:
            AssertionError if `get_nlaypg_max` does not match results
            calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        veg_trait_table = {
            1: {'nlaypg': 4},
            2: {'nlaypg': 6},
            3: {'nlaypg': 2},
        }

        known_nlaypg_max = 6
        nlaypg_max = forage.get_nlaypg_max(veg_trait_table)
        self.assertEqual(nlaypg_max, known_nlaypg_max)

    def test_distribute_water_to_soil_layer(self):
        """Test `distribute_water_to_soil_layer`.

        Use the function `distribute_water_to_soil_layer` to revise moisture
        content in one soil layer and calculate the moisture added to the next
        adjacent soil layer.

        Raises:
            AssertionError if `distribute_water_to_soil_layer` does not match
            value calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        array_size = (10, 10)
        tolerance = 0.00001

        # high moisture inputs, overflow to next soil layer
        adep = 13.68
        afiel = 0.32
        asmos = 3.1
        current_moisture_inputs = 8.291

        known_asmos_revised = 4.3776
        known_modified_moisture_inputs = 11.391

        adep_ar = numpy.full(array_size, adep)
        afiel_ar = numpy.full(array_size, afiel)
        asmos_ar = numpy.full(array_size, asmos)
        current_moisture_inputs_ar = numpy.full(
            array_size, current_moisture_inputs)

        insert_nodata_values_into_array(adep_ar, _IC_NODATA)
        insert_nodata_values_into_array(asmos_ar, _TARGET_NODATA)

        asmos_revised = forage.distribute_water_to_soil_layer(
            'asmos_revised')(
            adep_ar, afiel_ar, asmos_ar, current_moisture_inputs_ar)
        amov = forage.distribute_water_to_soil_layer(
            'amov')(adep_ar, afiel_ar, asmos_ar, current_moisture_inputs_ar)

        self.assert_all_values_in_array_within_range(
            asmos_revised, known_asmos_revised - tolerance,
            known_asmos_revised + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            amov, known_modified_moisture_inputs - tolerance,
            known_modified_moisture_inputs + tolerance, _TARGET_NODATA)

        # high field capacity, no overflow
        adep = 17.
        afiel = 0.482
        asmos = 0.01
        current_moisture_inputs = 4.2

        known_asmos_revised = 4.21
        known_modified_moisture_inputs = 0

        adep_ar = numpy.full(array_size, adep)
        afiel_ar = numpy.full(array_size, afiel)
        asmos_ar = numpy.full(array_size, asmos)
        current_moisture_inputs_ar = numpy.full(
            array_size, current_moisture_inputs)

        insert_nodata_values_into_array(afiel_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(
            current_moisture_inputs_ar, _TARGET_NODATA)

        asmos_revised = forage.distribute_water_to_soil_layer(
            'asmos_revised')(
            adep_ar, afiel_ar, asmos_ar, current_moisture_inputs_ar)
        amov = forage.distribute_water_to_soil_layer(
            'amov')(adep_ar, afiel_ar, asmos_ar, current_moisture_inputs_ar)

        self.assert_all_values_in_array_within_range(
            asmos_revised, known_asmos_revised - tolerance,
            known_asmos_revised + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            amov, known_modified_moisture_inputs - tolerance,
            known_modified_moisture_inputs + tolerance, _TARGET_NODATA)

    def test_calc_available_water_for_transpiration(self):
        """Test `calc_available_water_for_transpiration`.

        Use the function `calc_available_water_for_transpiration` to calculate
        water available in one soil layer for transpiration.

        Raises:
            AssertionError if `calc_available_water_for_transpiration` does not
            match values calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        array_size = (10, 10)
        tolerance = 0.0000001

        # low inputs, no water available for transpiration
        asmos = 3.6
        awilt = 0.52
        adep = 15
        known_avw = 0.

        asmos_ar = numpy.full(array_size, asmos)
        awilt_ar = numpy.full(array_size, awilt)
        adep_ar = numpy.full(array_size, adep)
        avw = forage.calc_available_water_for_transpiration(
            asmos_ar, awilt_ar, adep_ar)
        self.assert_all_values_in_array_within_range(
            avw, known_avw - tolerance, known_avw + tolerance, _TARGET_NODATA)

        # high moisture inputs, water available for transpiration
        asmos = 6.21
        awilt = 0.31
        adep = 15
        known_avw = 1.56

        asmos_ar = numpy.full(array_size, asmos)
        awilt_ar = numpy.full(array_size, awilt)
        adep_ar = numpy.full(array_size, adep)
        avw = forage.calc_available_water_for_transpiration(
            asmos_ar, awilt_ar, adep_ar)
        self.assert_all_values_in_array_within_range(
            avw, known_avw - tolerance, known_avw + tolerance, _TARGET_NODATA)

    def test_remove_transpiration(self):
        """ Test `remove_transpiration`.

        Use the function `remove_transpiration` to calculate asmos, revised
        moisture content of one soil layer, and avinj, water available for
        plant growth in the soil layer.

        Raises:
            AssertionError if `remove_transpiration` does not match value
                calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        array_size = (10, 10)
        tolerance = 0.00001

        # transpiration limited by current available moisture
        asmos = 4.15
        awilt = 0.26
        adep = 13.5
        trap = 4.17
        awwt = 1.428
        tot2 = 5.883

        known_asmos_revised = 3.51
        known_avinj = 0

        asmos_ar = numpy.full(array_size, asmos)
        awilt_ar = numpy.full(array_size, awilt)
        adep_ar = numpy.full(array_size, adep)
        trap_ar = numpy.full(array_size, trap)
        awwt_ar = numpy.full(array_size, awwt)
        tot2_ar = numpy.full(array_size, tot2)

        avinj = forage.remove_transpiration(
            'avinj')(asmos_ar, awilt_ar, adep_ar, trap_ar, awwt_ar, tot2_ar)
        asmos_revised = forage.remove_transpiration(
            'asmos')(asmos_ar, awilt_ar, adep_ar, trap_ar, awwt_ar, tot2_ar)
        self.assert_all_values_in_array_within_range(
            avinj, known_avinj - tolerance, known_avinj + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            asmos_revised, known_asmos_revised - tolerance,
            known_asmos_revised + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(asmos_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(awilt_ar, _TARGET_NODATA)

        avinj = forage.remove_transpiration(
            'avinj')(asmos_ar, awilt_ar, adep_ar, trap_ar, awwt_ar, tot2_ar)
        asmos_revised = forage.remove_transpiration(
            'asmos')(asmos_ar, awilt_ar, adep_ar, trap_ar, awwt_ar, tot2_ar)
        self.assert_all_values_in_array_within_range(
            avinj, known_avinj - tolerance, known_avinj + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            asmos_revised, known_asmos_revised - tolerance,
            known_asmos_revised + tolerance, _TARGET_NODATA)

        # transpiration limited by total potential transpiration
        asmos = 3.15
        awilt = 0.15
        adep = 12.5
        trap = 3.72
        awwt = 0.428
        tot2 = 4.883

        known_asmos_revised = 2.823938
        known_avinj = 0.948948

        asmos_ar = numpy.full(array_size, asmos)
        awilt_ar = numpy.full(array_size, awilt)
        adep_ar = numpy.full(array_size, adep)
        trap_ar = numpy.full(array_size, trap)
        awwt_ar = numpy.full(array_size, awwt)
        tot2_ar = numpy.full(array_size, tot2)

        avinj = forage.remove_transpiration(
            'avinj')(asmos_ar, awilt_ar, adep_ar, trap_ar, awwt_ar, tot2_ar)
        asmos_revised = forage.remove_transpiration(
            'asmos')(asmos_ar, awilt_ar, adep_ar, trap_ar, awwt_ar, tot2_ar)
        self.assert_all_values_in_array_within_range(
            avinj, known_avinj - tolerance, known_avinj + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            asmos_revised, known_asmos_revised - tolerance,
            known_asmos_revised + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(trap_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(adep_ar, _IC_NODATA)
        insert_nodata_values_into_array(tot2_ar, _TARGET_NODATA)

        avinj = forage.remove_transpiration(
            'avinj')(asmos_ar, awilt_ar, adep_ar, trap_ar, awwt_ar, tot2_ar)
        asmos_revised = forage.remove_transpiration(
            'asmos')(asmos_ar, awilt_ar, adep_ar, trap_ar, awwt_ar, tot2_ar)
        self.assert_all_values_in_array_within_range(
            avinj, known_avinj - tolerance, known_avinj + tolerance,
            _TARGET_NODATA)
        self.assert_all_values_in_array_within_range(
            asmos_revised, known_asmos_revised - tolerance,
            known_asmos_revised + tolerance, _TARGET_NODATA)

    def test_calc_relative_water_content_lyr_1(self):
        """Test `calc_relative_water_content_lyr_1`.

        Use the function `calc_relative_water_content_lyr_1` to calculate the
        relative water content of soil layer 1 prior to evaporative losses.

        Raises:
            AssertionError if `calc_relative_water_content_lyr_1` does not
                match values calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        array_size = (10, 10)
        tolerance = 0.00001

        asmos_1 = 2.52
        adep_1 = 16.42
        awilt_1 = 0.145
        afiel_1 = 0.774

        known_rwcf_1 = 0.013468

        asmos_1_ar = numpy.full(array_size, asmos_1)
        adep_1_ar = numpy.full(array_size, adep_1)
        awilt_1_ar = numpy.full(array_size, awilt_1)
        afiel_1_ar = numpy.full(array_size, afiel_1)

        rwcf_1 = forage.calc_relative_water_content_lyr_1(
            asmos_1_ar, adep_1_ar, awilt_1_ar, afiel_1_ar)
        self.assert_all_values_in_array_within_range(
            rwcf_1, known_rwcf_1 - tolerance, known_rwcf_1 + tolerance,
            _TARGET_NODATA)

        insert_nodata_values_into_array(asmos_1_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(afiel_1_ar, _TARGET_NODATA)

        rwcf_1 = forage.calc_relative_water_content_lyr_1(
            asmos_1_ar, adep_1_ar, awilt_1_ar, afiel_1_ar)
        self.assert_all_values_in_array_within_range(
            rwcf_1, known_rwcf_1 - tolerance, known_rwcf_1 + tolerance,
            _TARGET_NODATA)

    def test_calc_evaporation_loss(self):
        """Test `calc_evaporation_loss.

        Use the function `calc_evaporation_loss` to calculate moisture
        that evaporates from soil layer 1.

        Raises:
            AssertionError if `calc_evaporation_loss` does not match results
            calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        array_size = (10, 10)
        tolerance = 0.00001

        # limited by potential evaporation
        rwcf_1 = 0.99
        pevp = 3.1
        absevap = 2.1
        asmos_1 = 4.62
        awilt_1 = 0.153
        adep_1 = 14.2

        known_evlos = 0.64232

        rwcf_1_ar = numpy.full(array_size, rwcf_1)
        pevp_ar = numpy.full(array_size, pevp)
        absevap_ar = numpy.full(array_size, absevap)
        asmos_1_ar = numpy.full(array_size, asmos_1)
        awilt_1_ar = numpy.full(array_size, awilt_1)
        adep_1_ar = numpy.full(array_size, adep_1)

        evlos = forage.calc_evaporation_loss(
            rwcf_1_ar, pevp_ar, absevap_ar, asmos_1_ar, awilt_1_ar, adep_1_ar)
        self.assert_all_values_in_array_within_range(
            evlos, known_evlos - tolerance, known_evlos + tolerance,
            _TARGET_NODATA)

        insert_nodata_values_into_array(rwcf_1_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(asmos_1_ar, _SV_NODATA)

        evlos = forage.calc_evaporation_loss(
            rwcf_1_ar, pevp_ar, absevap_ar, asmos_1_ar, awilt_1_ar, adep_1_ar)
        self.assert_all_values_in_array_within_range(
            evlos, known_evlos - tolerance, known_evlos + tolerance,
            _TARGET_NODATA)

        # limited by available moisture
        rwcf_1 = 0.99
        pevp = 8.2
        absevap = 2.1
        asmos_1 = 3.6
        awilt_1 = 0.153
        adep_1 = 14.2

        known_evlos = 1.4274

        rwcf_1_ar = numpy.full(array_size, rwcf_1)
        pevp_ar = numpy.full(array_size, pevp)
        absevap_ar = numpy.full(array_size, absevap)
        asmos_1_ar = numpy.full(array_size, asmos_1)
        awilt_1_ar = numpy.full(array_size, awilt_1)
        adep_1_ar = numpy.full(array_size, adep_1)

        evlos = forage.calc_evaporation_loss(
            rwcf_1_ar, pevp_ar, absevap_ar, asmos_1_ar, awilt_1_ar, adep_1_ar)
        self.assert_all_values_in_array_within_range(
            evlos, known_evlos - tolerance, known_evlos + tolerance,
            _TARGET_NODATA)

        insert_nodata_values_into_array(pevp_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(adep_1_ar, _IC_NODATA)

        evlos = forage.calc_evaporation_loss(
            rwcf_1_ar, pevp_ar, absevap_ar, asmos_1_ar, awilt_1_ar, adep_1_ar)
        self.assert_all_values_in_array_within_range(
            evlos, known_evlos - tolerance, known_evlos + tolerance,
            _TARGET_NODATA)

    def test_raster_difference(self):
        """Test `raster_difference`.

        Use the function `raster_difference` to subtract one raster from
        another, while allowing nodata values in one raster to propagate to
        the result.  Then calculate the difference again, treating nodata
        values in the two rasters as zero.

        Raises:
            AssertionError if `raster_difference` does not match values
                calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage

        raster1_val = 10
        raster2_val = 3
        known_result = 7
        raster1_path = os.path.join(self.workspace_dir, 'raster1.tif')
        raster2_path = os.path.join(self.workspace_dir, 'raster2.tif')
        target_path = os.path.join(self.workspace_dir, 'target.tif')
        create_random_raster(raster1_path, raster1_val, raster1_val)
        create_random_raster(raster2_path, raster2_val, raster2_val)

        raster1_nodata = -99
        raster2_nodata = -999

        forage.raster_difference(
            raster1_path, raster1_nodata, raster2_path, raster2_nodata,
            target_path, _TARGET_NODATA, nodata_remove=True)
        self.assert_all_values_in_raster_within_range(
            target_path, known_result, known_result, _TARGET_NODATA)

        # rasters contain nodata, which should be propagated to result
        insert_nodata_values_into_raster(raster1_path, raster1_nodata)
        insert_nodata_values_into_raster(raster2_path, raster2_nodata)
        forage.raster_difference(
            raster1_path, raster1_nodata, raster2_path, raster2_nodata,
            target_path, _TARGET_NODATA, nodata_remove=False)
        self.assert_all_values_in_raster_within_range(
            target_path, known_result, known_result, _TARGET_NODATA)

        # full raster of nodata, which should be treated as zero
        create_random_raster(raster1_path, raster1_val, raster1_val)
        create_random_raster(raster2_path, raster2_nodata, raster2_nodata)
        forage.raster_difference(
            raster1_path, raster1_nodata, raster2_path, raster2_nodata,
            target_path, _TARGET_NODATA, nodata_remove=True)
        self.assert_all_values_in_raster_within_range(
            target_path, raster1_val, raster1_val, _TARGET_NODATA)

    def test_raster_sum(self):
        """Test `raster_sum`.

        Use the function `raster_sum` to add two rasters, while allowing
        nodata values in one raster to propagate to the result.  Then
        calculate the difference again, treating nodata values as zero.

        Raises:
            AssertionError if `raster_sum` does not match values calculated by
                hand

        Returns:
            None
        """
        from natcap.invest import forage

        raster1_val = 10
        raster2_val = 3
        known_result = raster1_val + raster2_val
        raster1_path = os.path.join(self.workspace_dir, 'raster1.tif')
        raster2_path = os.path.join(self.workspace_dir, 'raster2.tif')
        target_path = os.path.join(self.workspace_dir, 'target.tif')
        create_random_raster(raster1_path, raster1_val, raster1_val)
        create_random_raster(raster2_path, raster2_val, raster2_val)

        raster1_nodata = -99
        raster2_nodata = -999

        forage.raster_sum(
            raster1_path, raster1_nodata, raster2_path, raster2_nodata,
            target_path, _TARGET_NODATA, nodata_remove=True)
        self.assert_all_values_in_raster_within_range(
            target_path, known_result, known_result, _TARGET_NODATA)

        # rasters contain nodata, which should be propagated to result
        insert_nodata_values_into_raster(raster1_path, raster1_nodata)
        insert_nodata_values_into_raster(raster2_path, raster2_nodata)
        forage.raster_sum(
            raster1_path, raster1_nodata, raster2_path, raster2_nodata,
            target_path, _TARGET_NODATA, nodata_remove=False)
        self.assert_all_values_in_raster_within_range(
            target_path, known_result, known_result, _TARGET_NODATA)

        # full raster of nodata, which should be treated as zero
        create_random_raster(raster1_path, raster1_val, raster1_val)
        create_random_raster(raster2_path, raster2_nodata, raster2_nodata)
        forage.raster_sum(
            raster1_path, raster1_nodata, raster2_path, raster2_nodata,
            target_path, _TARGET_NODATA, nodata_remove=True)
        self.assert_all_values_in_raster_within_range(
            target_path, raster1_val, raster1_val, _TARGET_NODATA)

    def test_soil_water(self):
        """Test `soil_water`.

        Use the function `soil_water` to distribute precipitation inputs to
        snow, evaporation from snow, surface evaporation, and transpiration
        by plants. Compare results to values calculated by a point-based
        version defined here.

        Raises:
            AssertionError if `_soil_water` does not match values calculated
                by point-based version

        Returns:
            None
        """
        def soil_water_point(
                precip, max_temp, min_temp, snow, snlq, pet, tmelt_1, tmelt_2,
                shwave, strucc_1, metabc_1, fracro, precro, fwloss_1, fwloss_2,
                pft_dict, adep, afiel, awilt, awtl):
            """Point-based implementation of `soil_water`.

            Parameters:
                precip (float): precipitation for this month, cm
                max_temp (float): maximum average temperature for this month,
                    deg C
                min_temp (float): minimum average temperature for this month,
                    deg C
                snow (float): existing snowpack before new snowfall for this
                    month
                snlq (float): existing liquid in snowpack before new snowfall
                pet (float): reference evapotranspiration
                tmelt_1 (float): parameter, minimum temperature above which
                    snow will melt
                tmelt_2 (float): parameter, ratio between degrees above
                    the minimum temperature and cm of snow that will melt
                shwave (float): shortwave radiation outside the atmosphere
                strucc_1 (float): carbon in surface structural litter
                metabc_1 (float): metabolic surface carbon
                fracro (float):  parameter, fraction of surface water
                    above precro that is lost to runoff
                precro (float):  parameter, amount of surface water that
                    must be available for runoff to occur
                fwloss_1 (float):  parameter, scaling factor for interception
                    and evaporation of precip by vegetation
                fwloss_2 (float):  parameter, scaling factor for bare soil
                    evaporation of precip
                pft_dict (dict): dictionary containing parameters and % cover
                    for plant functional types
                adep (float): parameter, depth of each soil layer in cm
                afiel (float): field capacity of each soil layer
                awilt (float): wilting point of each soil layer
                awtl (float): parameter, transpiration weighting factor for
                    each soil layer

            Returns:
                dictionary of values:
                    snow, current snowpack
                    snlq, current liquid in snow
                    amov_2, water moving from soil layer 2
                    asmos_<lyr>, current soil moisture, for lyr in 1:nlaypg_max
                    avh2o_1_<PFT>, water available for plant growth, for PFT
                        in pft_id_set
                    avh2o_3, available soil moisture in top two soil layers
            """
            # snow
            tave = (max_temp + min_temp) / 2.
            inputs_after_snow = precip
            if tave <= 0:
                snow = snow + precip
                inputs_after_snow = 0
            else:
                if snow > 0:
                    snlq = snlq + precip
                    inputs_after_snow = 0
            if snow > 0:
                snowtot = snow + snlq
                evsnow = min(snowtot, (pet * 0.87))
                snow = max(snow - evsnow * (snow/snowtot), 0.)
                snlq = max(snlq-evsnow * (snlq/snowtot), 0.)
                pet_rem = max((pet - evsnow / 0.87), 0)

                if tave >= tmelt_1:
                    snowmelt = tmelt_2 * (tave - tmelt_1) * shwave
                    snowmelt = max(snowmelt, 0)
                    snowmelt = min(snowmelt, snow)
                    snow = snow - snowmelt
                    snlq = snlq + snowmelt
                    if snlq > (0.5 * snow):
                        add = snlq - 0.5 * snow
                        snlq = snlq - add
                        inputs_after_snow = add
            else:
                pet_rem = pet

            # canopy and litter cover influencing surface losses
            sum_aglivc = sum(
                [pft_dict[pft_i]['aglivc'] * pft_dict[pft_i]['cover'] for
                    pft_i in pft_dict.iterkeys()])
            sum_stdedc = sum(
                [pft_dict[pft_i]['stdedc'] * pft_dict[pft_i]['cover'] for
                    pft_i in pft_dict.iterkeys()])
            sum_tgprod = sum(
                [pft_dict[pft_i]['tgprod'] * pft_dict[pft_i]['cover'] for
                    pft_i in pft_dict.iterkeys()])
            alit = min((strucc_1 + metabc_1) * 2.5, 400)
            aliv = sum_aglivc * 2.5 + (0.25 * sum_tgprod)
            sd = min(aliv + (sum_stdedc * 2.5), 800.)
            # surface losses
            runoff = max(fracro * (inputs_after_snow - precro), 0.)
            inputs_after_surface = inputs_after_snow - runoff
            if snow == 0:
                aint = (0.0003 * alit + 0.0006 * sd) * fwloss_1
                absevap = (
                    0.5 * math.exp((-0.002 * alit) - (0.004 * sd)) * fwloss_2)
                evap_losses = min(
                    ((absevap + aint) * inputs_after_surface), 0.4 * pet_rem)
            else:
                absevap = 0
                evap_losses = 0
            inputs_after_surface = inputs_after_surface - evap_losses
            # potential transpiration
            current_moisture_inputs = inputs_after_surface
            trap = pet_rem - evap_losses
            if tave < 2:
                trap = 0
            else:
                trap = max(
                    min(trap, pet_rem * 0.65 *
                        (1 - math.exp(-0.02 * aliv))), 0)
            pevp = max(pet_rem - trap - evap_losses, 0.)
            tran = min(trap - 0.01, current_moisture_inputs)
            trap = trap - tran
            modified_moisture_inputs = current_moisture_inputs - tran
            # distribute moisture to soil layers prior to transpiration
            current_moisture_inputs = modified_moisture_inputs
            nlaypg_max = max(
                pft_dict[pft_i]['nlaypg'] for pft_i in pft_dict.iterkeys())
            asmos_dict = {lyr: asmos for lyr in xrange(1, nlaypg_max + 1)}
            amov_dict = {}
            for lyr in xrange(1, nlaypg_max + 1):
                afl = adep * afiel
                asmos_dict[lyr] = asmos_dict[lyr] + current_moisture_inputs
                if asmos_dict[lyr] > afl:
                    amov_dict[lyr] = asmos_dict[lyr]
                    asmos_dict[lyr] = afl
                else:
                    amov_dict[lyr] = 0
                current_moisture_inputs = amov_dict[lyr]
            # avw: available water for transpiration
            # awwt: water available for transpiration weighted by transpiration
            # depth for that soil layer
            awwt_dict = {}
            tot = 0
            tot2 = 0
            for lyr in xrange(1, nlaypg_max + 1):
                avw = max(asmos_dict[lyr] - awilt * adep, 0)
                awwt_dict[lyr] = avw * awtl
                tot = tot + avw
                tot2 = tot2 + awwt_dict[lyr]
            # revise total potential transpiration
            trap = min(trap, tot)
            # remove water via transpiration
            avinj_dict = {}
            for lyr in xrange(1, nlaypg_max + 1):
                avinj_dict[lyr] = max(asmos_dict[lyr] - awilt * adep, 0)
                trl = min((trap * awwt_dict[lyr]) / tot2, avinj_dict[lyr])
                avinj_dict[lyr] = avinj_dict[lyr] - trl
                asmos_dict[lyr] = asmos_dict[lyr] - trl
            # relative water content of soil layer 1
            rwcf_1 = (asmos_dict[1] / adep - awilt) / (afiel - awilt)
            # evaporation from soil layer 1
            evmt = max((rwcf_1 - 0.25) / (1 - 0.25), 0.01)
            evlos = min(
                evmt * pevp * absevap * 0.1,
                max(asmos_dict[1] - awilt * adep, 0))
            # remove evaporation from total moisture in soil layer 1
            asmos_dict[1] = asmos_dict[1] - evlos
            # remove evaporation from moisture available to plants in soil
            # layer 1
            avinj_dict[1] - avinj_dict[1] - evlos
            # calculate avh2o_1, soil water available for growth, for each PFT
            avh2o_1_dict = {}
            for pft_i in pft_dict.iterkeys():
                avh2o_1_dict[pft_i] = (
                    sum(avinj_dict[lyr] for lyr in
                        xrange(1, pft_dict[pft_i]['nlaypg'] + 1)) *
                    pft_dict[pft_i]['cover'])
            avh2o_3 = sum([avinj_dict[lyr] for lyr in [1, 2]])

            results_dict = {
                'snow': snow,
                'snlq': snlq,
                'amov_2': amov_dict[2],
                'avh2o_3': avh2o_3,
                'asmos': asmos_dict,
                'avh2o_1': avh2o_1_dict,
            }
            return results_dict

        def generate_model_inputs_from_point_inputs(
                precip, max_temp, min_temp, snow, snlq, pet, tmelt_1, tmelt_2,
                shwave, strucc_1, metabc_1, fracro, precro, fwloss_1, fwloss_2,
                pft_dict, adep, afiel, awilt, awtl):
            """Generate model inputs for `soil_water` from point inputs."""
            nrows = 1
            ncols = 1
            # aligned inputs
            aligned_inputs = {
                'max_temp_{}'.format(current_month): os.path.join(
                    self.workspace_dir,
                    'max_temp_{}.tif'.format(current_month)),
                'min_temp_{}'.format(current_month): os.path.join(
                    self.workspace_dir,
                    'min_temp_{}.tif'.format(current_month)),
                'precip_{}'.format(month_index): os.path.join(
                    self.workspace_dir, 'precip.tif'),
                'site_index': os.path.join(
                    self.workspace_dir, 'site_index.tif'),
            }
            create_random_raster(
                aligned_inputs['max_temp_{}'.format(current_month)], max_temp,
                max_temp, nrows=nrows, ncols=ncols)
            create_random_raster(
                aligned_inputs['min_temp_{}'.format(current_month)], min_temp,
                min_temp, nrows=nrows, ncols=ncols)
            create_random_raster(
                aligned_inputs['precip_{}'.format(month_index)], precip,
                precip, nrows=nrows, ncols=ncols)
            create_random_raster(
                aligned_inputs['site_index'], 1, 1, nrows=nrows, ncols=ncols)
            for pft_i in pft_dict.iterkeys():
                cover = pft_dict[pft_i]['cover']
                aligned_inputs['pft_{}'.format(pft_i)] = os.path.join(
                    self.workspace_dir, 'pft_{}.tif'.format(pft_i))
                create_random_raster(
                    aligned_inputs['pft_{}'.format(pft_i)], cover, cover,
                    nrows=nrows, ncols=ncols)

            # site param table
            site_param_table = {
                1: {
                    'tmelt_1': tmelt_1,
                    'tmelt_2': tmelt_2,
                    'fwloss_4': fwloss_4,
                    'fracro': fracro,
                    'precro': precro,
                    'fwloss_1': fwloss_1,
                    'fwloss_2': fwloss_2,
                }
            }
            for lyr in xrange(1, nlaypg_max + 1):
                site_param_table[1]['adep_{}'.format(lyr)] = adep
                site_param_table[1]['awtl_{}'.format(lyr)] = awtl
            # veg trait table
            pft_id_set = set([i for i in pft_dict.iterkeys()])
            veg_trait_table = {}
            for pft_i in pft_dict.iterkeys():
                veg_trait_table[pft_i] = {
                    'nlaypg': pft_dict[pft_i]['nlaypg'],
                }
            # previous state variables
            prev_sv_reg = {
                'strucc_1_path': os.path.join(
                    self.workspace_dir, 'strucc_1_prev.tif'),
                'metabc_1_path': os.path.join(
                    self.workspace_dir, 'metabc_1_prev.tif'),
                'snow_path': os.path.join(self.workspace_dir, 'snow_prev.tif'),
                'snlq_path': os.path.join(self.workspace_dir, 'snlq_prev.tif'),
            }
            create_random_raster(
                prev_sv_reg['strucc_1_path'], strucc_1, strucc_1, nrows=nrows,
                ncols=ncols)
            create_random_raster(
                prev_sv_reg['metabc_1_path'], metabc_1, metabc_1, nrows=nrows,
                ncols=ncols)
            create_random_raster(
                prev_sv_reg['snow_path'], snow, snow, nrows=nrows,
                ncols=ncols)
            create_random_raster(
                prev_sv_reg['snlq_path'], snlq, snlq, nrows=nrows,
                ncols=ncols)
            for lyr in xrange(1, nlaypg_max + 1):
                prev_sv_reg['asmos_{}_path'.format(lyr)] = os.path.join(
                    self.workspace_dir, 'asmos_{}_prev.tif'.format(lyr))
                create_random_raster(
                    prev_sv_reg['asmos_{}_path'.format(lyr)], asmos, asmos,
                    nrows=nrows, ncols=ncols)
            for pft_i in pft_dict.iterkeys():
                prev_sv_reg['aglivc_{}_path'.format(pft_i)] = os.path.join(
                    self.workspace_dir, 'aglivc_{}_prev.tif'.format(pft_i))
                prev_sv_reg['stdedc_{}_path'.format(pft_i)] = os.path.join(
                    self.workspace_dir, 'stdedc_{}_prev.tif'.format(pft_i))
                create_random_raster(
                    prev_sv_reg['aglivc_{}_path'.format(pft_i)],
                    pft_dict[pft_i]['aglivc'], pft_dict[pft_i]['aglivc'],
                    nrows=nrows, ncols=ncols)
                create_random_raster(
                    prev_sv_reg['stdedc_{}_path'.format(pft_i)],
                    pft_dict[pft_i]['stdedc'], pft_dict[pft_i]['stdedc'],
                    nrows=nrows, ncols=ncols)
            # current state variables
            sv_reg = {
                'snow_path': os.path.join(self.workspace_dir, 'snow.tif'),
                'snlq_path': os.path.join(self.workspace_dir, 'snlq.tif'),
                'avh2o_3_path': os.path.join(
                    self.workspace_dir, 'avh2o_3.tif'),
            }
            for lyr in xrange(1, nlaypg_max + 1):
                sv_reg['asmos_{}_path'.format(lyr)] = os.path.join(
                    self.workspace_dir, 'asmos_{}_path'.format(lyr))
            for pft_i in pft_dict.iterkeys():
                sv_reg['avh2o_1_{}_path'.format(pft_i)] = os.path.join(
                    self.workspace_dir, 'avh2o_1_{}.tif'.format(pft_i))
            # persistent parameters
            pp_reg = {}
            for lyr in xrange(1, nlaypg_max + 1):
                pp_reg['afiel_{}_path'.format(lyr)] = os.path.join(
                    self.workspace_dir, 'afiel_{}.tif'.format(lyr))
                pp_reg['awilt_{}_path'.format(lyr)] = os.path.join(
                    self.workspace_dir, 'awilt_{}.tif'.format(lyr))
                create_random_raster(
                    pp_reg['afiel_{}_path'.format(lyr)], afiel, afiel,
                    nrows=nrows, ncols=ncols)
                create_random_raster(
                    pp_reg['awilt_{}_path'.format(lyr)], awilt, awilt,
                    nrows=nrows, ncols=ncols)
            # monthly shared quantities
            month_reg = {
                'amov_2': os.path.join(self.workspace_dir, 'amov_2.tif'),
                'snowmelt': os.path.join(self.workspace_dir, 'snowmelt.tif')
            }
            for pft_i in pft_dict.iterkeys():
                month_reg['tgprod_{}'.format(pft_i)] = os.path.join(
                    self.workspace_dir, 'tgprod_{}.tif'.format(pft_i))
                create_random_raster(
                    month_reg['tgprod_{}'.format(pft_i)],
                    pft_dict[pft_i]['tgprod'], pft_dict[pft_i]['tgprod'],
                    nrows=nrows, ncols=ncols)
            input_dict = {
                'aligned_inputs': aligned_inputs,
                'site_param_table': site_param_table,
                'veg_trait_table': veg_trait_table,
                'current_month': current_month,
                'month_index': month_index,
                'prev_sv_reg': prev_sv_reg,
                'sv_reg': sv_reg,
                'pp_reg': pp_reg,
                'month_reg': month_reg,
                'pft_id_set': pft_id_set,
            }
            return input_dict
        from natcap.invest import forage

        # no snow, no snowfall
        pet = 4.9680004
        shwave = 457.95056

        # point inputs
        current_month = 10
        month_index = 9
        max_temp = 23.
        min_temp = -2.
        precip = 10.9

        # site parameters
        tmelt_1 = 0.
        tmelt_2 = 0.002
        fwloss_4 = 0.6
        fracro = 0.15
        precro = 8
        fwloss_1 = 0.8
        fwloss_2 = 0.779
        adep = 15
        awtl = 0.5

        # previous state variables
        strucc_1 = 46.1
        metabc_1 = 33.1
        snow = 0.
        snlq = 0.
        asmos = 2.04

        pft_dict = {
            1: {
                'cover': 0.4,
                'aglivc': 35.9,
                'stdedc': 49.874,
                'nlaypg': 5,
                'tgprod': 371.,
            },
            2: {
                'cover': 0.2,
                'aglivc': 52.1,
                'stdedc': 31.4,
                'nlaypg': 3,
                'tgprod': 300.2,
            },
            3: {
                'cover': 0.3,
                'aglivc': 20.77,
                'stdedc': 17.03,
                'nlaypg': 6,
                'tgprod': 200.8,
            },
        }
        nlaypg_max = max(
            pft_dict[pft_i]['nlaypg'] for pft_i in pft_dict.iterkeys())

        # persistent parameters
        afiel = 0.67
        awilt = 0.168

        tolerance = 0.000001
        results_dict = soil_water_point(
            precip, max_temp, min_temp, snow, snlq, pet, tmelt_1, tmelt_2,
            shwave, strucc_1, metabc_1, fracro, precro, fwloss_1, fwloss_2,
            pft_dict, adep, afiel, awilt, awtl)

        input_dict = generate_model_inputs_from_point_inputs(
            precip, max_temp, min_temp, snow, snlq, pet, tmelt_1, tmelt_2,
            shwave, strucc_1, metabc_1, fracro, precro, fwloss_1, fwloss_2,
            pft_dict, adep, afiel, awilt, awtl)
        forage._soil_water(
            input_dict['aligned_inputs'], input_dict['site_param_table'],
            input_dict['veg_trait_table'], input_dict['current_month'],
            input_dict['month_index'], input_dict['prev_sv_reg'],
            input_dict['sv_reg'], input_dict['pp_reg'],
            input_dict['month_reg'], input_dict['pft_id_set'])

        self.assert_all_values_in_raster_within_range(
            input_dict['sv_reg']['snow_path'],
            results_dict['snow'] - tolerance,
            results_dict['snow'] + tolerance, _SV_NODATA)
        self.assert_all_values_in_raster_within_range(
            input_dict['sv_reg']['snlq_path'],
            results_dict['snlq'] - tolerance,
            results_dict['snlq'] + tolerance, _SV_NODATA)
        for lyr in xrange(1, nlaypg_max + 1):
            self.assert_all_values_in_raster_within_range(
                input_dict['sv_reg']['asmos_{}_path'.format(lyr)],
                results_dict['asmos'][lyr] - tolerance,
                results_dict['asmos'][lyr] + tolerance, _SV_NODATA)
        self.assert_all_values_in_raster_within_range(
            input_dict['month_reg']['amov_2'],
            results_dict['amov_2'] - tolerance,
            results_dict['amov_2'] + tolerance, _TARGET_NODATA)
        for pft_i in input_dict['pft_id_set']:
            self.assert_all_values_in_raster_within_range(
                input_dict['sv_reg']['avh2o_1_{}_path'.format(pft_i)],
                results_dict['avh2o_1'][pft_i] - tolerance,
                results_dict['avh2o_1'][pft_i] + tolerance, _SV_NODATA)
        self.assert_all_values_in_raster_within_range(
            input_dict['sv_reg']['avh2o_3_path'],
            results_dict['avh2o_3'] - tolerance,
            results_dict['avh2o_3'] + tolerance, _SV_NODATA)

        # large snowmelt, large precip
        pet = 4.9680004
        shwave = 457.95056

        # point inputs
        current_month = 10
        month_index = 9
        max_temp = 23.
        min_temp = -2.
        precip = 30.117

        # site parameters
        tmelt_1 = 0.
        tmelt_2 = 0.002
        fwloss_4 = 0.6
        fracro = 0.15
        precro = 8
        fwloss_1 = 0.8
        fwloss_2 = 0.779
        adep = 15
        awtl = 0.5

        # previous state variables
        strucc_1 = 46.1
        metabc_1 = 33.1
        snow = 7.2
        snlq = snow / 2.
        asmos = 0.882

        pft_dict = {
            1: {
                'cover': 0.4,
                'aglivc': 85.9,
                'stdedc': 49.874,
                'nlaypg': 5,
                'tgprod': 371.,
            },
            2: {
                'cover': 0.2,
                'aglivc': 12.1,
                'stdedc': 31.4,
                'nlaypg': 3,
                'tgprod': 300.2,
            },
            3: {
                'cover': 0.3,
                'aglivc': 27.77,
                'stdedc': 17.03,
                'nlaypg': 6,
                'tgprod': 200.8,
            },
        }
        nlaypg_max = max(
            pft_dict[pft_i]['nlaypg'] for pft_i in pft_dict.iterkeys())

        # persistent parameters
        afiel = 0.67
        awilt = 0.168

        tolerance = 0.00001
        amov_tolerance = 0.011
        results_dict = soil_water_point(
            precip, max_temp, min_temp, snow, snlq, pet, tmelt_1, tmelt_2,
            shwave, strucc_1, metabc_1, fracro, precro, fwloss_1, fwloss_2,
            pft_dict, adep, afiel, awilt, awtl)

        input_dict = generate_model_inputs_from_point_inputs(
            precip, max_temp, min_temp, snow, snlq, pet, tmelt_1, tmelt_2,
            shwave, strucc_1, metabc_1, fracro, precro, fwloss_1, fwloss_2,
            pft_dict, adep, afiel, awilt, awtl)
        forage._soil_water(
            input_dict['aligned_inputs'], input_dict['site_param_table'],
            input_dict['veg_trait_table'], input_dict['current_month'],
            input_dict['month_index'], input_dict['prev_sv_reg'],
            input_dict['sv_reg'], input_dict['pp_reg'],
            input_dict['month_reg'], input_dict['pft_id_set'])

        self.assert_all_values_in_raster_within_range(
            input_dict['sv_reg']['snow_path'],
            results_dict['snow'] - tolerance,
            results_dict['snow'] + tolerance, _SV_NODATA)
        self.assert_all_values_in_raster_within_range(
            input_dict['sv_reg']['snlq_path'],
            results_dict['snlq'] - tolerance,
            results_dict['snlq'] + tolerance, _SV_NODATA)
        for lyr in xrange(1, nlaypg_max + 1):
            self.assert_all_values_in_raster_within_range(
                input_dict['sv_reg']['asmos_{}_path'.format(lyr)],
                results_dict['asmos'][lyr] - tolerance,
                results_dict['asmos'][lyr] + tolerance, _SV_NODATA)
        self.assert_all_values_in_raster_within_range(
            input_dict['month_reg']['amov_2'],
            results_dict['amov_2'] - amov_tolerance,
            results_dict['amov_2'] + amov_tolerance, _TARGET_NODATA)
        for pft_i in input_dict['pft_id_set']:
            self.assert_all_values_in_raster_within_range(
                input_dict['sv_reg']['avh2o_1_{}_path'.format(pft_i)],
                results_dict['avh2o_1'][pft_i] - tolerance,
                results_dict['avh2o_1'][pft_i] + tolerance, _TARGET_NODATA)
        self.assert_all_values_in_raster_within_range(
            input_dict['sv_reg']['avh2o_3_path'],
            results_dict['avh2o_3'] - tolerance,
            results_dict['avh2o_3'] + tolerance, _SV_NODATA)

    def test_monthly_N_fixation(self):
        """Test `_monthly_N_fixation`.

        Use the function `_monthly_N_fixation` to calculate monthly atmospheric
        N deposition and add it to the surface mineral N pool.  Compare results
        to values calculated by point-based version.

        Raises:
            AssertionError if updated mineral_1_1 does not match value
                calculated by point-based version

        Returns:
            None
        """
        from natcap.invest import forage

        # known inputs
        precip = 12.6
        month_index = 4
        annual_precip = 230.
        baseNdep = 24.
        epnfs_2 = 0.01
        prev_minerl_1_1 = 3.2
        minerl_1_1 = monthly_N_fixation_point(
            precip, annual_precip, baseNdep, epnfs_2, prev_minerl_1_1)

        aligned_inputs = {
            'precip_{}'.format(month_index): os.path.join(
                self.workspace_dir, 'precip_{}.tif'.format(month_index)),
            'site_index': os.path.join(
                    self.workspace_dir, 'site_index.tif'),
        }
        year_reg = {
            'annual_precip_path': os.path.join(
                self.workspace_dir, 'annual_precip.tif'),
            'baseNdep_path': os.path.join(self.workspace_dir, 'baseNdep.tif')
        }
        prev_sv_reg = {
            'minerl_1_1_path': os.path.join(
                self.workspace_dir, 'minerl_1_1_prev.tif'),
        }
        sv_reg = {
            'minerl_1_1_path': os.path.join(
                self.workspace_dir, 'minerl_1_1.tif'),
        }
        site_param_table = {1: {'epnfs_2': epnfs_2}}

        create_random_raster(
            aligned_inputs['precip_{}'.format(month_index)], precip,
            precip)
        create_random_raster(aligned_inputs['site_index'], 1, 1)
        create_random_raster(
            year_reg['annual_precip_path'], annual_precip, annual_precip)
        create_random_raster(year_reg['baseNdep_path'], baseNdep, baseNdep)
        create_random_raster(
            prev_sv_reg['minerl_1_1_path'], prev_minerl_1_1, prev_minerl_1_1)

        tolerance = 0.000001
        forage._monthly_N_fixation(
            aligned_inputs, month_index, site_param_table,
            year_reg, prev_sv_reg, sv_reg)
        self.assert_all_values_in_raster_within_range(
            sv_reg['minerl_1_1_path'], minerl_1_1 - tolerance,
            minerl_1_1 + tolerance, _SV_NODATA)

        insert_nodata_values_into_raster(
            aligned_inputs['precip_{}'.format(month_index)], -9999)
        insert_nodata_values_into_raster(
            year_reg['baseNdep_path'], _TARGET_NODATA)
        insert_nodata_values_into_raster(
            prev_sv_reg['minerl_1_1_path'], _SV_NODATA)

        forage._monthly_N_fixation(
            aligned_inputs, month_index, site_param_table,
            year_reg, prev_sv_reg, sv_reg)
        self.assert_all_values_in_raster_within_range(
            sv_reg['minerl_1_1_path'], minerl_1_1 - tolerance,
            minerl_1_1 + tolerance, _SV_NODATA)

    def test_calc_anerb(self):
        """Test `calc_anerb`.

        Use the function `calc_anerb` to calculate the effect of soil anaerobic
        conditions on decomposition. Compare the calculated value to a value
        calculated by point-based version.

        Raises:
            AssertionError if anerb does not match value calculated by
                point-based version

        Returns:
            None
        """

        from natcap.invest import forage

        array_shape = (10, 10)
        tolerance = 0.00000001

        # low rprpet, anerb = 1
        rprpet = 0.8824
        pevap = 6.061683
        drain = 0.003
        aneref_1 = 1.5
        aneref_2 = 3.
        aneref_3 = 0.3

        rprpet_arr = numpy.full(array_shape, rprpet)
        pevap_arr = numpy.full(array_shape, pevap)
        drain_arr = numpy.full(array_shape, drain)
        aneref_1_arr = numpy.full(array_shape, aneref_1)
        aneref_2_arr = numpy.full(array_shape, aneref_2)
        aneref_3_arr = numpy.full(array_shape, aneref_3)

        anerb = calc_anerb_point(
            rprpet, pevap, drain, aneref_1, aneref_2, aneref_3)
        anerb_arr = forage.calc_anerb(
            rprpet_arr, pevap_arr, drain_arr, aneref_1_arr, aneref_2_arr,
            aneref_3_arr)
        self.assert_all_values_in_array_within_range(
            anerb_arr, anerb - tolerance, anerb + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(rprpet_arr, _TARGET_NODATA)
        anerb_arr = forage.calc_anerb(
            rprpet_arr, pevap_arr, drain_arr, aneref_1_arr, aneref_2_arr,
            aneref_3_arr)
        self.assert_all_values_in_array_within_range(
            anerb_arr, anerb - tolerance, anerb + tolerance, _TARGET_NODATA)

        # high rprpet, xh2o > 0
        rprpet = 2.0004
        rprpet_arr = numpy.full(array_shape, rprpet)
        anerb = calc_anerb_point(
            rprpet, pevap, drain, aneref_1, aneref_2, aneref_3)
        anerb_arr = forage.calc_anerb(
            rprpet_arr, pevap_arr, drain_arr, aneref_1_arr, aneref_2_arr,
            aneref_3_arr)
        self.assert_all_values_in_array_within_range(
            anerb_arr, anerb - tolerance, anerb + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(drain_arr, _IC_NODATA)
        anerb_arr = forage.calc_anerb(
            rprpet_arr, pevap_arr, drain_arr, aneref_1_arr, aneref_2_arr,
            aneref_3_arr)
        self.assert_all_values_in_array_within_range(
            anerb_arr, anerb - tolerance, anerb + tolerance, _TARGET_NODATA)

        # high rprpet, xh2o = 0
        drain = 1.
        drain_arr = numpy.full(array_shape, drain)
        anerb = calc_anerb_point(
            rprpet, pevap, drain, aneref_1, aneref_2, aneref_3)
        anerb_arr = forage.calc_anerb(
            rprpet_arr, pevap_arr, drain_arr, aneref_1_arr, aneref_2_arr,
            aneref_3_arr)
        self.assert_all_values_in_array_within_range(
            anerb_arr, anerb - tolerance, anerb + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(pevap_arr, _TARGET_NODATA)
        anerb_arr = forage.calc_anerb(
            rprpet_arr, pevap_arr, drain_arr, aneref_1_arr, aneref_2_arr,
            aneref_3_arr)
        self.assert_all_values_in_array_within_range(
            anerb_arr, anerb - tolerance, anerb + tolerance, _TARGET_NODATA)

    def test_declig_point(self):
        """Test `declig_point`.

        Use the function `declig_point` to calculate change in state variables
        as material containing lignin decomposes into SOM2 and SOM1. Compare
        calculated changes in state variables to changes calculated by hand.

        Raises:
            AssertionError if `declig_point` does not match values calculated
                by hand

        Returns:
            None
        """
        # no decomposition happens, mineral ratios are insufficient
        aminrl_1 = 0.
        aminrl_2 = 0.
        ligcon = 0.3779
        rsplig = 0.0146
        ps1co2_lyr = 0.0883
        strucc_lyr = 155.5253
        tcflow = 40.82
        struce_lyr_1 = 0.7776
        struce_lyr_2 = 0.3111
        rnew_lyr_1_1 = 190.  # observed ratio: 200.0068
        rnew_lyr_2_1 = 300.  # observed ratio: 499.9206
        rnew_lyr_1_2 = 0.
        rnew_lyr_2_2 = 0.
        minerl_1_1 = 0.
        minerl_1_2 = 0.

        d_strucc_lyr = declig_point(
            'd_strucc')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_struce_lyr_1 = declig_point(
            'd_struce_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_struce_lyr_2 = declig_point(
            'd_struce_2')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_minerl_1_1 = declig_point(
            'd_minerl_1_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_minerl_1_2 = declig_point(
            'd_minerl_1_2')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_gromin_1 = declig_point(
            'd_gromin_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som2c_lyr = declig_point(
            'd_som2c')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som2e_lyr_1 = declig_point(
            'd_som2e_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som2e_lyr_2 = declig_point(
            'd_som2e_2')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som1c_lyr = declig_point(
            'd_som1c')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som1e_lyr_1 = declig_point(
            'd_som1e_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som1e_lyr_2 = declig_point(
            'd_som1e_2')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)

        self.assertEqual(d_strucc_lyr, 0)
        self.assertEqual(d_struce_lyr_1, 0)
        self.assertEqual(d_struce_lyr_2, 0)
        self.assertEqual(d_minerl_1_1, 0)
        self.assertEqual(d_minerl_1_2, 0)
        self.assertEqual(d_gromin_1, 0)
        self.assertEqual(d_som2c_lyr, 0)
        self.assertEqual(d_som2e_lyr_1, 0)
        self.assertEqual(d_som2e_lyr_2, 0)
        self.assertEqual(d_som1c_lyr, 0)
        self.assertEqual(d_som1e_lyr_1, 0)
        self.assertEqual(d_som1e_lyr_2, 0)

        # decomposition occurs, subsidized by mineral N and P
        aminrl_1 = 6.4944
        aminrl_2 = 33.2791
        ligcon = 0.3779
        rsplig = 0.0146
        ps1co2_lyr = 0.0883
        strucc_lyr = 155.5253
        tcflow = 40.82
        struce_lyr_1 = 0.7776
        struce_lyr_2 = 0.3111
        rnew_lyr_1_1 = 210.8  # greater than observed ratio: 200.0068
        rnew_lyr_2_1 = 540.2  # greater than observed ratio: 499.9206
        rnew_lyr_1_2 = 190.3
        rnew_lyr_2_2 = 520.8
        minerl_1_1 = 6.01
        minerl_1_2 = 32.87

        # values calculated by hand
        d_strucc_lyr_obs = -40.82
        d_struce_lyr_1_obs = -0.204093044668617
        d_struce_lyr_2_obs = -0.08165296578756
        d_minerl_1_1_obs = 0.014387319170996
        d_minerl_1_2_obs = 0.00960796090459662
        d_gromin_1_obs = 0.0182639608121622
        d_som2c_lyr_obs = 15.2006601812
        d_som2e_lyr_1_obs = 0.0798773525023647
        d_som2e_lyr_2_obs = 0.029187135524578
        d_som1c_lyr_obs = 23.1518210274
        d_som1e_lyr_1_obs = 0.109828372995256
        d_som1e_lyr_2_obs = 0.0428578693583858

        # call 12 times, 12 return values
        d_strucc_lyr = declig_point(
            'd_strucc')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_struce_lyr_1 = declig_point(
            'd_struce_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_struce_lyr_2 = declig_point(
            'd_struce_2')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_minerl_1_1 = declig_point(
            'd_minerl_1_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_minerl_1_2 = declig_point(
            'd_minerl_1_2')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_gromin_1 = declig_point(
            'd_gromin_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som2c_lyr = declig_point(
            'd_som2c')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som2e_lyr_1 = declig_point(
            'd_som2e_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som2e_lyr_2 = declig_point(
            'd_som2e_2')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som1c_lyr = declig_point(
            'd_som1c')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som1e_lyr_1 = declig_point(
            'd_som1e_1')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)
        d_som1e_lyr_2 = declig_point(
            'd_som1e_2')(
                aminrl_1, aminrl_2, ligcon, rsplig, ps1co2_lyr, strucc_lyr,
                tcflow, struce_lyr_1, struce_lyr_2, rnew_lyr_1_1, rnew_lyr_2_1,
                rnew_lyr_1_2, rnew_lyr_2_2, minerl_1_1, minerl_1_2)

        self.assertAlmostEqual(d_strucc_lyr, d_strucc_lyr_obs, places=10)
        self.assertAlmostEqual(d_struce_lyr_1, d_struce_lyr_1_obs, places=10)
        self.assertAlmostEqual(d_struce_lyr_2, d_struce_lyr_2_obs, places=10)
        self.assertAlmostEqual(d_minerl_1_1, d_minerl_1_1_obs, places=10)
        self.assertAlmostEqual(d_minerl_1_2, d_minerl_1_2_obs, places=10)
        self.assertAlmostEqual(d_gromin_1, d_gromin_1_obs, places=10)
        self.assertAlmostEqual(d_som2c_lyr, d_som2c_lyr_obs, places=10)
        self.assertAlmostEqual(d_som2e_lyr_1, d_som2e_lyr_1_obs, places=10)
        self.assertAlmostEqual(d_som2e_lyr_2, d_som2e_lyr_2_obs, places=10)
        self.assertAlmostEqual(d_som1c_lyr, d_som1c_lyr_obs, places=10)
        self.assertAlmostEqual(d_som1e_lyr_1, d_som1e_lyr_1_obs, places=10)
        self.assertAlmostEqual(d_som1e_lyr_2, d_som1e_lyr_2_obs, places=10)

    def test_decomposition(self):
        """Test `_decomposition`.

        Test the function `_decomposition` to calculate the change in state
        variables following one decomposition step.  Compare modified state
        variables to values calculated by point version of decomposition.

        Raises:
            AssertionError if change in state variable calculated by
                `_decomposition` does not match change calculated by point-
                based version

        Returns:
            None
        """
        def generate_model_inputs_from_point_inputs(
                inputs, params, state_var_dict, year_reg_vals, month_reg_vals,
                pp_reg_vals, rnew_dict):
            """Generate model inputs for `_decomposition` from test inputs."""
            nrows = 1
            ncols = 1

            aligned_inputs = {
                'precip_{}'.format(month_index): os.path.join(
                    self.workspace_dir, 'precip.tif'),
                'min_temp_{}'.format(current_month): os.path.join(
                    self.workspace_dir, 'min_temp.tif'),
                'max_temp_{}'.format(current_month): os.path.join(
                    self.workspace_dir, 'max_temp.tif'),
                'ph_path': os.path.join(self.workspace_dir, 'pH.tif'),
                'site_index': os.path.join(
                    self.workspace_dir, 'site_index.tif'),
            }
            create_random_raster(
                aligned_inputs['precip_{}'.format(month_index)],
                inputs['precip'], inputs['precip'], nrows=nrows, ncols=ncols)
            create_random_raster(
                aligned_inputs['min_temp_{}'.format(current_month)],
                inputs['min_temp'], inputs['min_temp'],
                nrows=nrows, ncols=ncols)
            create_random_raster(
                aligned_inputs['max_temp_{}'.format(current_month)],
                inputs['max_temp'], inputs['max_temp'], nrows=nrows,
                ncols=ncols)
            create_random_raster(
                aligned_inputs['ph_path'], inputs['pH'], inputs['pH'],
                nrows=nrows, ncols=ncols)
            create_random_raster(
                aligned_inputs['site_index'], 1, 1, nrows=nrows, ncols=ncols)

            site_param_table = {1: {}}
            for key, value in params.iteritems():
                site_param_table[1][key] = value

            year_reg = {
                'annual_precip_path': os.path.join(
                    self.workspace_dir, 'annual_precip.tif'),
                'baseNdep_path': os.path.join(
                    self.workspace_dir, 'baseNdep.tif'),
            }
            create_random_raster(
                year_reg['annual_precip_path'], year_reg_vals['annual_precip'],
                year_reg_vals['annual_precip'], nrows=nrows, ncols=ncols)
            create_random_raster(
                year_reg['baseNdep_path'], year_reg_vals['baseNdep'],
                year_reg_vals['baseNdep'], nrows=nrows, ncols=ncols)

            month_reg = {
                'snowmelt': os.path.join(self.workspace_dir, 'snowmelt.tif'),
                'amov_2': os.path.join(self.workspace_dir, 'amov_2.tif'),
            }
            create_random_raster(
                month_reg['snowmelt'], month_reg_vals['snowmelt'],
                month_reg_vals['snowmelt'], nrows=nrows, ncols=ncols)
            create_random_raster(
                month_reg['amov_2'], month_reg_vals['amov_2'],
                month_reg_vals['amov_2'], nrows=nrows, ncols=ncols)

            prev_sv_reg = {}
            sv_reg = {}
            for state_var in [
                    'minerl_1_1', 'minerl_1_2', 'snow', 'avh2o_3', 'parent_2',
                    'secndy_2', 'occlud']:
                prev_sv_reg['{}_path'.format(state_var)] = os.path.join(
                    self.workspace_dir, '{}_p.tif'.format(state_var))
                sv_reg['{}_path'.format(state_var)] = os.path.join(
                    self.workspace_dir, '{}.tif'.format(state_var))
                create_random_raster(
                    prev_sv_reg['{}_path'.format(state_var)],
                    state_var_dict[state_var], state_var_dict[state_var],
                    nrows=nrows, ncols=ncols)
                create_random_raster(
                    sv_reg['{}_path'.format(state_var)],
                    state_var_dict[state_var], state_var_dict[state_var],
                    nrows=nrows, ncols=ncols)
            for compartment in ['strlig']:
                for lyr in [1, 2]:
                    state_var = '{}_{}'.format(compartment, lyr)
                    prev_sv_reg['{}_path'.format(state_var)] = os.path.join(
                        self.workspace_dir, '{}_p.tif'.format(state_var))
                    sv_reg['{}_path'.format(state_var)] = os.path.join(
                        self.workspace_dir, '{}.tif'.format(state_var))
                    create_random_raster(
                        prev_sv_reg['{}_path'.format(state_var)],
                        state_var_dict[state_var], state_var_dict[state_var],
                        nrows=nrows, ncols=ncols)
            for compartment in ['som3']:
                state_var = '{}c'.format(compartment)
                prev_sv_reg['{}_path'.format(state_var)] = os.path.join(
                    self.workspace_dir, '{}_p.tif'.format(state_var))
                sv_reg['{}_path'.format(state_var)] = os.path.join(
                    self.workspace_dir, '{}.tif'.format(state_var))
                create_random_raster(
                    prev_sv_reg['{}_path'.format(state_var)],
                    state_var_dict[state_var],
                    state_var_dict[state_var],
                    nrows=nrows, ncols=ncols)
                for iel in [1, 2]:
                    state_var = '{}e_{}'.format(compartment, iel)
                    prev_sv_reg['{}_path'.format(state_var)] = os.path.join(
                        self.workspace_dir, '{}_p.tif'.format(state_var))
                    sv_reg['{}_path'.format(state_var)] = os.path.join(
                        self.workspace_dir, '{}.tif'.format(state_var))
                    create_random_raster(
                        prev_sv_reg['{}_path'.format(state_var)],
                        state_var_dict[state_var],
                        state_var_dict[state_var],
                        nrows=nrows, ncols=ncols)
            for compartment in ['struc', 'metab', 'som1', 'som2']:
                for lyr in [1, 2]:
                    state_var = '{}c_{}'.format(compartment, lyr)
                    prev_sv_reg['{}_path'.format(state_var)] = os.path.join(
                        self.workspace_dir, '{}_p.tif'.format(state_var))
                    sv_reg['{}_path'.format(state_var)] = os.path.join(
                        self.workspace_dir, '{}.tif'.format(state_var))
                    create_random_raster(
                        prev_sv_reg['{}_path'.format(state_var)],
                        state_var_dict[state_var],
                        state_var_dict[state_var],
                        nrows=nrows, ncols=ncols)
                    for iel in [1, 2]:
                        state_var = '{}e_{}_{}'.format(compartment, lyr, iel)
                        prev_sv_reg['{}_path'.format(state_var)] = (
                            os.path.join(
                                self.workspace_dir,
                                '{}_p.tif'.format(state_var)))
                        sv_reg['{}_path'.format(state_var)] = os.path.join(
                            self.workspace_dir, '{}.tif'.format(state_var))
                        create_random_raster(
                            prev_sv_reg['{}_path'.format(state_var)],
                            state_var_dict[state_var],
                            state_var_dict[state_var],
                            nrows=nrows, ncols=ncols)
            pp_reg = {
                'rnewas_1_1': os.path.join(
                    self.workspace_dir, 'rnewas_1_1.tif'),
                'rnewas_1_2': os.path.join(
                    self.workspace_dir, 'rnewas_1_2.tif'),
                'rnewas_2_1': os.path.join(
                    self.workspace_dir, 'rnewas_2_1.tif'),
                'rnewas_2_2': os.path.join(
                    self.workspace_dir, 'rnewas_2_2.tif'),
                'rnewbs_1_1': os.path.join(
                    self.workspace_dir, 'rnewbs_1_1.tif'),
                'rnewbs_1_2': os.path.join(
                    self.workspace_dir, 'rnewbs_1_2.tif'),
                'rnewbs_2_1': os.path.join(
                    self.workspace_dir, 'rnewbs_2_1.tif'),
                'rnewbs_2_2': os.path.join(
                    self.workspace_dir, 'rnewbs_2_2.tif'),
            }
            for key in rnew_dict.iterkeys():
                create_random_raster(
                    pp_reg[key], rnew_dict[key], rnew_dict[key],
                    nrows=nrows, ncols=ncols)
            pp_reg['eftext_path'] = os.path.join(
                self.workspace_dir, 'eftext.tif')
            pp_reg['orglch_path'] = os.path.join(
                self.workspace_dir, 'orglch.tif')
            pp_reg['fps1s3_path'] = os.path.join(
                self.workspace_dir, 'fps1s3.tif')
            pp_reg['fps2s3_path'] = os.path.join(
                self.workspace_dir, 'fps2s3.tif')
            create_random_raster(
                pp_reg['eftext_path'], pp_reg_vals['eftext'],
                pp_reg_vals['eftext'], nrows=nrows, ncols=ncols)
            create_random_raster(
                pp_reg['orglch_path'], pp_reg_vals['orglch'],
                pp_reg_vals['orglch'], nrows=nrows, ncols=ncols)
            create_random_raster(
                pp_reg['fps1s3_path'], pp_reg_vals['fps1s3'],
                pp_reg_vals['fps1s3'], nrows=nrows, ncols=ncols)
            create_random_raster(
                pp_reg['fps2s3_path'], pp_reg_vals['fps2s3'],
                pp_reg_vals['fps2s3'], nrows=nrows, ncols=ncols)

            input_dict = {
                'aligned_inputs': aligned_inputs,
                'site_param_table': site_param_table,
                'year_reg': year_reg,
                'month_reg': month_reg,
                'prev_sv_reg': prev_sv_reg,
                'sv_reg': sv_reg,
                'pp_reg': pp_reg,
            }
            return input_dict
        from natcap.invest import forage

        # known inputs
        current_month = 4
        month_index = 2
        inputs = {
            'min_temp': 3.2,
            'max_temp': 17.73,
            'precip': 30.5,
            'pH': 6.84,
        }
        params = {
            'fwloss_4': 0.8,
            'teff_1': 15.4,
            'teff_2': 11.75,
            'teff_3': 29.7,
            'teff_4': 0.031,
            'drain': 1,
            'aneref_1': 1.5,
            'aneref_2': 3,
            'aneref_3': 0.3,
            'epnfs_2': 30.,
            'sorpmx': 2,
            'pslsrb': 1,
            'strmax_1': 5000,
            'dec1_1': 3.9,
            'pligst_1': 3,
            'rsplig': 0.3,
            'ps1co2_1': 0.45,
            'strmax_2': 5000,
            'dec1_2': 4.9,
            'pligst_2': 3,
            'ps1co2_2': 0.55,
            'pcemic1_1_1': 16,
            'pcemic1_2_1': 10,
            'pcemic1_3_1': 0.02,
            'pcemic1_1_2': 200,
            'pcemic1_2_2': 99,
            'pcemic1_3_2': 0.0015,
            'dec2_1': 14.8,
            'pmco2_1': 0.55,
            'varat1_1_1': 14,
            'varat1_2_1': 3,
            'varat1_3_1': 2,
            'varat1_1_2': 150,
            'varat1_2_2': 30,
            'varat1_3_2': 2,
            'dec2_2': 18.5,
            'pmco2_2': 0.55,
            'rad1p_1_1': 12,
            'rad1p_2_1': 3,
            'rad1p_3_1': 5,
            'rad1p_1_2': 220,
            'rad1p_2_2': 5,
            'rad1p_3_2': 100,
            'dec3_1': 6,
            'p1co2a_1': 0.6,
            'varat22_1_1': 20,
            'varat22_2_1': 12,
            'varat22_3_1': 2,
            'varat22_1_2': 400,
            'varat22_2_2': 100,
            'varat22_3_2': 2,
            'dec3_2': 7.3,
            'p1co2_2': 0.55,
            'animpt': 5,
            'varat3_1_1': 8,
            'varat3_2_1': 6,
            'varat3_3_1': 2,
            'varat3_1_2': 200,
            'varat3_2_2': 50,
            'varat3_3_2': 2,
            'omlech_3': 60,
            'dec5_2': 0.2,
            'p2co2_2': 0.55,
            'dec5_1': 0.2,
            'p2co2_1': 0.55,
            'dec4': 0.0045,
            'p3co2': 0.55,
            'cmix': 0.5,
            'pparmn_2': 0.0001,
            'psecmn_2': 0.0022,
            'nlayer': 5,
            'pmnsec_2': 0,
            'psecoc1': 0,
            'psecoc2': 0,
            'vlossg': 1,
        }
        state_var_dict = {
            'snow': 0.,
            'avh2o_3': 3.110,
            'strlig_1': 0.3779,
            'strlig_2': 0.2871,
            'minerl_1_1': 6.4143,
            'minerl_1_2': 33.1954,
            'strucc_1': 156.0546,
            'struce_1_1': 0.7803,
            'struce_1_2': 0.3121,
            'som2c_1': 92.9935,
            'som2e_1_1': 4.2466,
            'som2e_1_2': 0.1328,
            'som1c_1': 12.8192,
            'som1e_1_1': 1.5752,
            'som1e_1_2': 0.1328,
            'strucc_2': 163.2008,
            'struce_2_1': 0.816,
            'struce_2_2': 0.3264,
            'som2c_2': 922.1382,
            'som2e_2_1': 60.0671,
            'som2e_2_2': 5.6575,
            'som1c_2': 39.1055,
            'som1e_2_1': 12.2767,
            'som1e_2_2': 1.2263,
            'metabc_1': 13.7579,
            'metabe_1_1': 1.1972,
            'metabe_1_2': 0.0351,
            'metabc_2': 9.8169,
            'metabe_2_1': 0.6309,
            'metabe_2_2': 0.088,
            'som3c': 544.8848,
            'som3e_1': 87.4508,
            'som3e_2': 79.917,
            'parent_2': 100.,
            'secndy_2': 50.,
            'occlud': 50.,
        }
        for lyr in xrange(1, params['nlayer'] + 1):
            state_var_dict['minerl_{}_2'.format(lyr)] = 20.29

        year_reg_vals = {
            'annual_precip': 230.,
            'baseNdep': 24.,
        }
        month_reg_vals = {
            'snowmelt': 0.29,
            'amov_2': 0.481,
        }
        pp_reg_vals = {
            'eftext': 0.15,
            'orglch': 0.01,
            'fps1s3': 0.58,
            'fps2s3': 0.58,
        }
        rnew_dict = {
            'rnewas_1_1': 210.8,
            'rnewas_2_1': 540.2,
            'rnewas_1_2': 190.3,
            'rnewas_2_2': 520.8,
            'rnewbs_1_1': 210.8,
            'rnewbs_2_1': 540.2,
            'rnewbs_1_2': 190.3,
            'rnewbs_2_2': 520.8,
        }
        pevap = 9.324202
        input_dict = generate_model_inputs_from_point_inputs(
            inputs, params, state_var_dict, year_reg_vals, month_reg_vals,
            pp_reg_vals, rnew_dict)

        sv_mod_point = decomposition_point(
            inputs, params, state_var_dict, year_reg_vals, month_reg_vals,
            pp_reg_vals, rnew_dict, pevap)

        sv_mod_raster = forage._decomposition(
            input_dict['aligned_inputs'], current_month, month_index,
            input_dict['site_param_table'], input_dict['year_reg'],
            input_dict['month_reg'], input_dict['prev_sv_reg'],
            input_dict['sv_reg'], input_dict['pp_reg'])

        for compartment in ['struc', 'metab', 'som1', 'som2']:
            tolerance = 0.0001
            for lyr in [1, 2]:
                state_var = '{}c_{}'.format(compartment, lyr)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
                for iel in [1, 2]:
                    state_var = '{}e_{}_{}'.format(compartment, lyr, iel)
                    point_value = sv_mod_point[state_var]
                    model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                    self.assert_all_values_in_raster_within_range(
                        model_res_path, point_value - tolerance,
                        point_value + tolerance, _SV_NODATA)
        for compartment in ['som3']:
            state_var = '{}c'.format(compartment)
            point_value = sv_mod_point[state_var]
            model_res_path = sv_mod_raster['{}_path'.format(state_var)]
            self.assert_all_values_in_raster_within_range(
                model_res_path, point_value - tolerance,
                point_value + tolerance, _SV_NODATA)
            for iel in [1, 2]:
                state_var = '{}e_{}'.format(compartment, iel)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
        for compartment in ['minerl_1']:
            tolerance = 0.001
            for iel in [1, 2]:
                state_var = '{}_{}'.format(compartment, iel)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
        for state_var in ['parent_2', 'secndy_2', 'occlud']:
            point_value = sv_mod_point[state_var]
            model_res_path = sv_mod_raster['{}_path'.format(state_var)]
            self.assert_all_values_in_raster_within_range(
                model_res_path, point_value - tolerance,
                point_value + tolerance, _SV_NODATA)

        # no decomposition, mineral ratios are insufficient
        rnew_dict = {
            'rnewas_1_1': 190.,
            'rnewas_2_1': 300.,
            'rnewas_1_2': 140.,
            'rnewas_2_2': 200.,
            'rnewbs_1_1': 190.,
            'rnewbs_2_1': 300.,
            'rnewbs_1_2': 140.,
            'rnewbs_2_2': 200.
        }
        state_var_dict['minerl_1_1'] = 0.00000001
        state_var_dict['minerl_1_2'] = 0.00000001

        input_dict = generate_model_inputs_from_point_inputs(
            inputs, params, state_var_dict, year_reg_vals, month_reg_vals,
            pp_reg_vals, rnew_dict)

        sv_mod_point = decomposition_point(
            inputs, params, state_var_dict, year_reg_vals, month_reg_vals,
            pp_reg_vals, rnew_dict, pevap)

        sv_mod_raster = forage._decomposition(
            input_dict['aligned_inputs'], current_month, month_index,
            input_dict['site_param_table'], input_dict['year_reg'],
            input_dict['month_reg'], input_dict['prev_sv_reg'],
            input_dict['sv_reg'], input_dict['pp_reg'])

        for compartment in ['struc', 'metab', 'som1', 'som2']:
            tolerance = 0.0001
            for lyr in [1, 2]:
                state_var = '{}c_{}'.format(compartment, lyr)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
                for iel in [1, 2]:
                    state_var = '{}e_{}_{}'.format(compartment, lyr, iel)
                    point_value = sv_mod_point[state_var]
                    model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                    self.assert_all_values_in_raster_within_range(
                        model_res_path, point_value - tolerance,
                        point_value + tolerance, _SV_NODATA)
        for compartment in ['som3']:
            state_var = '{}c'.format(compartment)
            point_value = sv_mod_point[state_var]
            model_res_path = sv_mod_raster['{}_path'.format(state_var)]
            self.assert_all_values_in_raster_within_range(
                model_res_path, point_value - tolerance,
                point_value + tolerance, _SV_NODATA)
            for iel in [1, 2]:
                state_var = '{}e_{}'.format(compartment, iel)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
        for compartment in ['minerl_1']:
            tolerance = 0.001
            for iel in [1, 2]:
                state_var = '{}_{}'.format(compartment, iel)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
        for state_var in ['parent_2', 'secndy_2', 'occlud']:
            point_value = sv_mod_point[state_var]
            model_res_path = sv_mod_raster['{}_path'.format(state_var)]
            self.assert_all_values_in_raster_within_range(
                model_res_path, point_value - tolerance,
                point_value + tolerance, _SV_NODATA)

        # decomposition occurs, subsidized by mineral N and P
        rnew_dict = {
            'rnewas_1_1': 210.8,
            'rnewas_2_1': 540.2,
            'rnewas_1_2': 190.3,
            'rnewas_2_2': 520.8,
            'rnewbs_1_1': 210.8,
            'rnewbs_2_1': 540.2,
            'rnewbs_1_2': 190.3,
            'rnewbs_2_2': 520.8
        }
        state_var_dict['minerl_1_1'] = 6.01
        state_var_dict['minerl_1_2'] = 32.87

        state_var_dict['strlig_1'] = 0.2987
        state_var_dict['strlig_2'] = 0.4992

        input_dict = generate_model_inputs_from_point_inputs(
            inputs, params, state_var_dict, year_reg_vals, month_reg_vals,
            pp_reg_vals, rnew_dict)

        sv_mod_point = decomposition_point(
            inputs, params, state_var_dict, year_reg_vals, month_reg_vals,
            pp_reg_vals, rnew_dict, pevap)

        sv_mod_raster = forage._decomposition(
            input_dict['aligned_inputs'], current_month, month_index,
            input_dict['site_param_table'], input_dict['year_reg'],
            input_dict['month_reg'], input_dict['prev_sv_reg'],
            input_dict['sv_reg'], input_dict['pp_reg'])

        for compartment in ['struc', 'metab', 'som1', 'som2']:
            tolerance = 0.0001
            for lyr in [1, 2]:
                state_var = '{}c_{}'.format(compartment, lyr)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
                for iel in [1, 2]:
                    state_var = '{}e_{}_{}'.format(compartment, lyr, iel)
                    point_value = sv_mod_point[state_var]
                    model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                    self.assert_all_values_in_raster_within_range(
                        model_res_path, point_value - tolerance,
                        point_value + tolerance, _SV_NODATA)
        for compartment in ['som3']:
            state_var = '{}c'.format(compartment)
            point_value = sv_mod_point[state_var]
            model_res_path = sv_mod_raster['{}_path'.format(state_var)]
            self.assert_all_values_in_raster_within_range(
                model_res_path, point_value - tolerance,
                point_value + tolerance, _SV_NODATA)
            for iel in [1, 2]:
                state_var = '{}e_{}'.format(compartment, iel)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
        for compartment in ['minerl_1']:
            tolerance = 0.001
            for iel in [1, 2]:
                state_var = '{}_{}'.format(compartment, iel)
                point_value = sv_mod_point[state_var]
                model_res_path = sv_mod_raster['{}_path'.format(state_var)]
                self.assert_all_values_in_raster_within_range(
                    model_res_path, point_value - tolerance,
                    point_value + tolerance, _SV_NODATA)
        for state_var in ['parent_2', 'secndy_2', 'occlud']:
            point_value = sv_mod_point[state_var]
            model_res_path = sv_mod_raster['{}_path'.format(state_var)]
            self.assert_all_values_in_raster_within_range(
                model_res_path, point_value - tolerance,
                point_value + tolerance, _SV_NODATA)

    def test_esched(self):
        """Test `esched`.

        Use the function `esched` to calculate the flow of one element
        accompanying decomposition of C.  Test that the function matches
        values calculated by point-based version.

        Raises:
            AssertionError if `esched` does not match value calculated
                by `esched_point`

        Returns:
            None
        """
        from natcap.invest import forage
        tolerance = 0.00000001

        # immobilization
        cflow = 15.2006
        tca = 155.5253
        rcetob = 190.3
        anps = 0.7776
        labile = 6.01

        material_leaving_a = esched_point(
            'material_leaving_a')(cflow, tca, rcetob, anps, labile)
        material_arriving_b = esched_point(
            'material_arriving_b')(cflow, tca, rcetob, anps, labile)
        mineral_flow = esched_point(
            'mineral_flow')(cflow, tca, rcetob, anps, labile)

        # raster inputs
        cflow_path = os.path.join(self.workspace_dir, 'cflow.tif')
        tca_path = os.path.join(self.workspace_dir, 'tca.tif')
        rcetob_path = os.path.join(self.workspace_dir, 'rcetob.tif')
        anps_path = os.path.join(self.workspace_dir, 'anps.tif')
        labile_path = os.path.join(self.workspace_dir, 'labile.tif')
        # output paths
        mat_leaving_a_path = os.path.join(self.workspace_dir, 'leavinga.tif')
        mat_arriving_b_path = os.path.join(self.workspace_dir, 'arrivingb.tif')
        mineral_flow_path = os.path.join(self.workspace_dir, 'mineralflow.tif')

        create_random_raster(cflow_path, cflow, cflow)
        create_random_raster(tca_path, tca, tca)
        create_random_raster(rcetob_path, rcetob, rcetob)
        create_random_raster(anps_path, anps, anps)
        create_random_raster(labile_path, labile, labile)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_leaving_a'), mat_leaving_a_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_arriving_b'), mat_arriving_b_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('mineral_flow'), mineral_flow_path,
            gdal.GDT_Float32, _IC_NODATA)

        self.assert_all_values_in_raster_within_range(
            mat_leaving_a_path, material_leaving_a - tolerance,
            material_leaving_a + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mat_arriving_b_path, material_arriving_b - tolerance,
            material_arriving_b + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mineral_flow_path, mineral_flow - tolerance,
            mineral_flow + tolerance, _IC_NODATA)

        insert_nodata_values_into_raster(cflow_path, _IC_NODATA)
        insert_nodata_values_into_raster(anps_path, _SV_NODATA)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_leaving_a'), mat_leaving_a_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_arriving_b'), mat_arriving_b_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('mineral_flow'), mineral_flow_path,
            gdal.GDT_Float32, _IC_NODATA)

        self.assert_all_values_in_raster_within_range(
            mat_leaving_a_path, material_leaving_a - tolerance,
            material_leaving_a + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mat_arriving_b_path, material_arriving_b - tolerance,
            material_arriving_b + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mineral_flow_path, mineral_flow - tolerance,
            mineral_flow + tolerance, _IC_NODATA)

        # mineralization
        cflow = 15.2006
        tca = 155.5253
        rcetob = 520.8
        anps = 0.3111
        labile = 32.87

        material_leaving_a = esched_point(
            'material_leaving_a')(cflow, tca, rcetob, anps, labile)
        material_arriving_b = esched_point(
            'material_arriving_b')(cflow, tca, rcetob, anps, labile)
        mineral_flow = esched_point(
            'mineral_flow')(cflow, tca, rcetob, anps, labile)

        create_random_raster(cflow_path, cflow, cflow)
        create_random_raster(tca_path, tca, tca)
        create_random_raster(rcetob_path, rcetob, rcetob)
        create_random_raster(anps_path, anps, anps)
        create_random_raster(labile_path, labile, labile)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_leaving_a'), mat_leaving_a_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_arriving_b'), mat_arriving_b_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('mineral_flow'), mineral_flow_path,
            gdal.GDT_Float32, _IC_NODATA)

        self.assert_all_values_in_raster_within_range(
            mat_leaving_a_path, material_leaving_a - tolerance,
            material_leaving_a + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mat_arriving_b_path, material_arriving_b - tolerance,
            material_arriving_b + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mineral_flow_path, mineral_flow - tolerance,
            mineral_flow + tolerance, _IC_NODATA)

        insert_nodata_values_into_raster(anps_path, _SV_NODATA)
        insert_nodata_values_into_raster(tca_path, _SV_NODATA)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_leaving_a'), mat_leaving_a_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_arriving_b'), mat_arriving_b_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('mineral_flow'), mineral_flow_path,
            gdal.GDT_Float32, _IC_NODATA)

        self.assert_all_values_in_raster_within_range(
            mat_leaving_a_path, material_leaving_a - tolerance,
            material_leaving_a + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mat_arriving_b_path, material_arriving_b - tolerance,
            material_arriving_b + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mineral_flow_path, mineral_flow - tolerance,
            mineral_flow + tolerance, _IC_NODATA)

        # no movement
        cflow = 15.2006
        tca = 155.5253
        rcetob = 520.8
        anps = 0.3111
        labile = 0.

        material_leaving_a = esched_point(
            'material_leaving_a')(cflow, tca, rcetob, anps, labile)
        material_arriving_b = esched_point(
            'material_arriving_b')(cflow, tca, rcetob, anps, labile)
        mineral_flow = esched_point(
            'mineral_flow')(cflow, tca, rcetob, anps, labile)

        create_random_raster(cflow_path, cflow, cflow)
        create_random_raster(tca_path, tca, tca)
        create_random_raster(rcetob_path, rcetob, rcetob)
        create_random_raster(anps_path, anps, anps)
        create_random_raster(labile_path, labile, labile)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_leaving_a'), mat_leaving_a_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_arriving_b'), mat_arriving_b_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('mineral_flow'), mineral_flow_path,
            gdal.GDT_Float32, _IC_NODATA)

        self.assert_all_values_in_raster_within_range(
            mat_leaving_a_path, material_leaving_a - tolerance,
            material_leaving_a + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mat_arriving_b_path, material_arriving_b - tolerance,
            material_arriving_b + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mineral_flow_path, mineral_flow - tolerance,
            mineral_flow + tolerance, _IC_NODATA)

        insert_nodata_values_into_raster(rcetob_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(labile_path, _SV_NODATA)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_leaving_a'), mat_leaving_a_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('material_arriving_b'), mat_arriving_b_path,
            gdal.GDT_Float32, _IC_NODATA)
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                cflow_path, tca_path, rcetob_path, anps_path, labile_path]],
            forage.esched('mineral_flow'), mineral_flow_path,
            gdal.GDT_Float32, _IC_NODATA)

        self.assert_all_values_in_raster_within_range(
            mat_leaving_a_path, material_leaving_a - tolerance,
            material_leaving_a + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mat_arriving_b_path, material_arriving_b - tolerance,
            material_arriving_b + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            mineral_flow_path, mineral_flow - tolerance,
            mineral_flow + tolerance, _IC_NODATA)

    def test_nutrient_flow(self):
        """Test `nutrient_flow`.

        Use the function `nutrient_flow` to calculate and apply the flow of
        one element accompanying decomposition of C. Test calculated values
        against values calculated by point-based version.

        Raises:
            AssertionError if `nutrient_flow` does not match values calculated
                by `esched_point`

        Returns:
            None
        """
        from natcap.invest import forage
        tolerance = 0.00000001

        # immobilization
        cflow = 15.2006
        tca = 155.5253
        rcetob = 190.3
        anps = 0.7776
        labile = 6.01

        material_leaving_a = esched_point(
            'material_leaving_a')(cflow, tca, rcetob, anps, labile)
        material_arriving_b = esched_point(
            'material_arriving_b')(cflow, tca, rcetob, anps, labile)
        mineral_flow = esched_point(
            'mineral_flow')(cflow, tca, rcetob, anps, labile)

        d_estatv_donating = -material_leaving_a
        d_estatv_receiving = material_arriving_b
        d_minerl = mineral_flow
        if mineral_flow > 0:
            gromin = mineral_flow
        else:
            gromin = 0

        # raster inputs
        cflow_path = os.path.join(self.workspace_dir, 'cflow.tif')
        tca_path = os.path.join(self.workspace_dir, 'tca.tif')
        rcetob_path = os.path.join(self.workspace_dir, 'rcetob.tif')
        anps_path = os.path.join(self.workspace_dir, 'anps.tif')
        labile_path = os.path.join(self.workspace_dir, 'labile.tif')
        d_estatv_donating_path = os.path.join(
            self.workspace_dir, 'estatv_donating.tif')
        d_estatv_receiving_path = os.path.join(
            self.workspace_dir, 'estatv_receiving.tif')
        d_minerl_path = os.path.join(self.workspace_dir, 'minerl.tif')
        gromin_path = os.path.join(self.workspace_dir, 'gromin.tif')

        create_random_raster(cflow_path, cflow, cflow)
        create_random_raster(tca_path, tca, tca)
        create_random_raster(rcetob_path, rcetob, rcetob)
        create_random_raster(anps_path, anps, anps)
        create_random_raster(labile_path, labile, labile)
        create_random_raster(d_estatv_donating_path, 0, 0)
        create_random_raster(d_estatv_receiving_path, 0, 0)
        create_random_raster(d_minerl_path, 0, 0)
        create_random_raster(gromin_path, 0, 0)

        forage.nutrient_flow(
            cflow_path, tca_path, anps_path, rcetob_path,
            labile_path, d_estatv_donating_path, d_estatv_receiving_path,
            d_minerl_path, gromin_path)
        self.assert_all_values_in_raster_within_range(
            d_estatv_donating_path, d_estatv_donating - tolerance,
            d_estatv_donating + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            d_estatv_receiving_path, d_estatv_receiving - tolerance,
            d_estatv_receiving + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            d_minerl_path, d_minerl - tolerance, d_minerl + tolerance,
            _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            gromin_path, gromin - tolerance, gromin + tolerance, _IC_NODATA)

        create_random_raster(d_estatv_donating_path, 0, 0)
        create_random_raster(d_estatv_receiving_path, 0, 0)
        create_random_raster(d_minerl_path, 0, 0)
        create_random_raster(gromin_path, 0, 0)

        insert_nodata_values_into_raster(cflow_path, _IC_NODATA)
        insert_nodata_values_into_raster(tca_path, _SV_NODATA)
        insert_nodata_values_into_raster(rcetob_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(anps_path, _SV_NODATA)
        insert_nodata_values_into_raster(d_estatv_donating_path, _IC_NODATA)

        forage.nutrient_flow(
            cflow_path, tca_path, anps_path, rcetob_path,
            labile_path, d_estatv_donating_path, d_estatv_receiving_path,
            d_minerl_path, gromin_path)
        self.assert_all_values_in_raster_within_range(
            d_estatv_donating_path, d_estatv_donating - tolerance,
            d_estatv_donating + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            d_estatv_receiving_path, d_estatv_receiving - tolerance,
            d_estatv_receiving + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            d_minerl_path, d_minerl - tolerance, d_minerl + tolerance,
            _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            gromin_path, gromin - tolerance, gromin + tolerance,
            _TARGET_NODATA)

        # mineralization
        cflow = 15.2006
        tca = 155.5253
        rcetob = 520.8
        anps = 0.3111
        labile = 32.87

        material_leaving_a = esched_point(
            'material_leaving_a')(cflow, tca, rcetob, anps, labile)
        material_arriving_b = esched_point(
            'material_arriving_b')(cflow, tca, rcetob, anps, labile)
        mineral_flow = esched_point(
            'mineral_flow')(cflow, tca, rcetob, anps, labile)

        d_estatv_donating = -material_leaving_a
        d_estatv_receiving = material_arriving_b
        d_minerl = mineral_flow

        # raster inputs
        create_random_raster(cflow_path, cflow, cflow)
        create_random_raster(tca_path, tca, tca)
        create_random_raster(rcetob_path, rcetob, rcetob)
        create_random_raster(anps_path, anps, anps)
        create_random_raster(labile_path, labile, labile)
        create_random_raster(d_estatv_donating_path, 0, 0)
        create_random_raster(d_estatv_receiving_path, 0, 0)
        create_random_raster(d_minerl_path, 0, 0)

        forage.nutrient_flow(
            cflow_path, tca_path, anps_path, rcetob_path,
            labile_path, d_estatv_donating_path, d_estatv_receiving_path,
            d_minerl_path)
        self.assert_all_values_in_raster_within_range(
            d_estatv_donating_path, d_estatv_donating - tolerance,
            d_estatv_donating + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            d_estatv_receiving_path, d_estatv_receiving - tolerance,
            d_estatv_receiving + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            d_minerl_path, d_minerl - tolerance, d_minerl + tolerance,
            _IC_NODATA)

        create_random_raster(d_estatv_donating_path, 0, 0)
        create_random_raster(d_estatv_receiving_path, 0, 0)
        create_random_raster(d_minerl_path, 0, 0)

        insert_nodata_values_into_raster(cflow_path, _IC_NODATA)
        insert_nodata_values_into_raster(tca_path, _SV_NODATA)
        insert_nodata_values_into_raster(rcetob_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(anps_path, _SV_NODATA)
        insert_nodata_values_into_raster(d_estatv_donating_path, _IC_NODATA)

        forage.nutrient_flow(
            cflow_path, tca_path, anps_path, rcetob_path,
            labile_path, d_estatv_donating_path, d_estatv_receiving_path,
            d_minerl_path)
        self.assert_all_values_in_raster_within_range(
            d_estatv_donating_path, d_estatv_donating - tolerance,
            d_estatv_donating + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            d_estatv_receiving_path, d_estatv_receiving - tolerance,
            d_estatv_receiving + tolerance, _IC_NODATA)
        self.assert_all_values_in_raster_within_range(
            d_minerl_path, d_minerl - tolerance, d_minerl + tolerance,
            _IC_NODATA)

    def test_fsfunc(self):
        """Test `fsfunc`.

        Use the function `fsfunc` to calculate the fraction of mineral P that
        is in solution.  Compare calculated value to value calculated by point
        version.

        Raises:
            AssertionError if `initialize_aminrl_2` does not match value
                calculated by point version of function

        Returns:
            None
        """
        from natcap.invest import forage

        tolerance = 0.00001

        # known values
        minerl_1_2 = 32.87
        sorpmx = 2.
        pslsrb = 1.
        fsol_point = fsfunc_point(minerl_1_2, pslsrb, sorpmx)

        # raster inputs
        minerl_1_2_path = os.path.join(self.workspace_dir, 'minerl_1_2.tif')
        sorpmx_path = os.path.join(self.workspace_dir, 'sorpmx.tif')
        pslsrb_path = os.path.join(self.workspace_dir, 'pslsrb.tif')
        fsol_path = os.path.join(self.workspace_dir, 'fsol.tif')

        create_random_raster(minerl_1_2_path, minerl_1_2, minerl_1_2)
        create_random_raster(sorpmx_path, sorpmx, sorpmx)
        create_random_raster(pslsrb_path, pslsrb, pslsrb)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                minerl_1_2_path, sorpmx_path, pslsrb_path]],
            forage.fsfunc, fsol_path, gdal.GDT_Float32,
            _SV_NODATA)
        self.assert_all_values_in_raster_within_range(
            fsol_path, fsol_point - tolerance, fsol_point + tolerance,
            _SV_NODATA)

        insert_nodata_values_into_raster(minerl_1_2_path, _SV_NODATA)
        insert_nodata_values_into_raster(pslsrb_path, _IC_NODATA)

        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                minerl_1_2_path, sorpmx_path, pslsrb_path]],
            forage.fsfunc, fsol_path, gdal.GDT_Float32,
            _SV_NODATA)
        self.assert_all_values_in_raster_within_range(
            fsol_path, fsol_point - tolerance, fsol_point + tolerance,
            _SV_NODATA)

    def test_calc_tcflow_strucc_1(self):
        """Test `calc_tcflow_strucc_1`.

        Use the function `calc_tcflow_strucc_1` to calculate total flow
        out of surface structural C. Ensure that calculated values match
        values calculated by point-based version defined here.

        Raises:
            AssertionError if `calc_tcflow_strucc_1` does not match values
                calculated by point-based version

        Returns:
            None
        """
        def tcflow_strucc_1_point(
                aminrl_1, aminrl_2, strucc_1, struce_1_1, struce_1_2,
                rnewas_1_1, rnewas_2_1, strmax_1, defac, dec1_1, pligst_1,
                strlig_1, pheff_struc):
            """Point-based implementation of `calc_tcflow_strucc_1`.

            Returns:
                tcflow_strucc_1, total flow of C limited by N and P
            """
            potential_flow = (min(
                strucc_1, strmax_1) * defac * dec1_1 *
                math.exp(-pligst_1 * strlig_1) * 0.020833 * pheff_struc)

            decompose_mask = (
                ((aminrl_1 > 0.0000001) | (
                    (strucc_1 / struce_1_1) <= rnewas_1_1)) &
                ((aminrl_2 > 0.0000001) | (
                    (strucc_1 / struce_1_2) <= rnewas_2_1)))

            if decompose_mask:
                tcflow_strucc_1 = potential_flow
            else:
                tcflow_strucc_1 = 0
            return tcflow_strucc_1
        from natcap.invest import forage

        array_shape = (10, 10)
        tolerance = 0.0000001

        # decomposition can occur
        aminrl_1 = 6.4143
        aminrl_2 = 30.9253
        strucc_1 = 156.0546
        struce_1_1 = 0.7803
        struce_1_2 = 0.3121
        rnewas_1_1 = 210.8
        rnewas_2_1 = 540.2
        strmax_1 = 5000.
        defac = 0.822
        dec1_1 = 3.9
        pligst_1 = 3.
        strlig_1 = 0.3779
        pH = 6.84
        pheff_struc = numpy.clip(
            (0.5 + (1.1 / numpy.pi) *
                numpy.arctan(numpy.pi * 0.7 * (pH - 4.))), 0, 1)

        tcflow_strucc_1 = tcflow_strucc_1_point(
            aminrl_1, aminrl_2, strucc_1, struce_1_1, struce_1_2,
            rnewas_1_1, rnewas_2_1, strmax_1, defac, dec1_1, pligst_1,
            strlig_1, pheff_struc)

        # array inputs
        aminrl_1_ar = numpy.full(array_shape, aminrl_1)
        aminrl_2_ar = numpy.full(array_shape, aminrl_2)
        strucc_1_ar = numpy.full(array_shape, strucc_1)
        struce_1_1_ar = numpy.full(array_shape, struce_1_1)
        struce_1_2_ar = numpy.full(array_shape, struce_1_2)
        rnewas_1_1_ar = numpy.full(array_shape, rnewas_1_1)
        rnewas_2_1_ar = numpy.full(array_shape, rnewas_2_1)
        strmax_1_ar = numpy.full(array_shape, strmax_1)
        defac_ar = numpy.full(array_shape, defac)
        dec1_1_ar = numpy.full(array_shape, dec1_1)
        pligst_1_ar = numpy.full(array_shape, pligst_1)
        strlig_1_ar = numpy.full(array_shape, strlig_1)
        pheff_struc_ar = numpy.full(array_shape, pheff_struc)

        tcflow_strucc1_ar = forage.calc_tcflow_strucc_1(
            aminrl_1_ar, aminrl_2_ar, strucc_1_ar, struce_1_1_ar,
            struce_1_2_ar, rnewas_1_1_ar, rnewas_2_1_ar, strmax_1_ar, defac_ar,
            dec1_1_ar, pligst_1_ar, strlig_1_ar, pheff_struc_ar)

        self.assert_all_values_in_array_within_range(
            tcflow_strucc1_ar, tcflow_strucc_1 - tolerance,
            tcflow_strucc_1 + tolerance, _IC_NODATA)

        insert_nodata_values_into_array(struce_1_2_ar, _SV_NODATA)
        insert_nodata_values_into_array(defac_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(strlig_1_ar, _SV_NODATA)

        tcflow_strucc1_ar = forage.calc_tcflow_strucc_1(
            aminrl_1_ar, aminrl_2_ar, strucc_1_ar, struce_1_1_ar,
            struce_1_2_ar, rnewas_1_1_ar, rnewas_2_1_ar, strmax_1_ar, defac_ar,
            dec1_1_ar, pligst_1_ar, strlig_1_ar, pheff_struc_ar)

        self.assert_all_values_in_array_within_range(
            tcflow_strucc1_ar, tcflow_strucc_1 - tolerance,
            tcflow_strucc_1 + tolerance, _IC_NODATA)

        # N insufficient to allow decomposition
        aminrl_1 = 0.
        aminrl_2 = 30.9253
        strucc_1 = 156.0546
        struce_1_1 = 0.7803
        struce_1_2 = 0.3121
        rnewas_1_1 = 170.
        rnewas_2_1 = 540.2
        strmax_1 = 5000.
        defac = 0.822
        dec1_1 = 3.9
        pligst_1 = 3.
        strlig_1 = 0.3779
        pH = 6.84
        pheff_struc = numpy.clip(
            (0.5 + (1.1 / numpy.pi) *
                numpy.arctan(numpy.pi * 0.7 * (pH - 4.))), 0, 1)

        tcflow_strucc_1 = tcflow_strucc_1_point(
            aminrl_1, aminrl_2, strucc_1, struce_1_1, struce_1_2,
            rnewas_1_1, rnewas_2_1, strmax_1, defac, dec1_1, pligst_1,
            strlig_1, pheff_struc)

        # array inputs
        aminrl_1_ar = numpy.full(array_shape, aminrl_1)
        aminrl_2_ar = numpy.full(array_shape, aminrl_2)
        strucc_1_ar = numpy.full(array_shape, strucc_1)
        struce_1_1_ar = numpy.full(array_shape, struce_1_1)
        struce_1_2_ar = numpy.full(array_shape, struce_1_2)
        rnewas_1_1_ar = numpy.full(array_shape, rnewas_1_1)
        rnewas_2_1_ar = numpy.full(array_shape, rnewas_2_1)
        strmax_1_ar = numpy.full(array_shape, strmax_1)
        defac_ar = numpy.full(array_shape, defac)
        dec1_1_ar = numpy.full(array_shape, dec1_1)
        pligst_1_ar = numpy.full(array_shape, pligst_1)
        strlig_1_ar = numpy.full(array_shape, strlig_1)
        pheff_struc_ar = numpy.full(array_shape, pheff_struc)

        tcflow_strucc1_ar = forage.calc_tcflow_strucc_1(
            aminrl_1_ar, aminrl_2_ar, strucc_1_ar, struce_1_1_ar,
            struce_1_2_ar, rnewas_1_1_ar, rnewas_2_1_ar, strmax_1_ar, defac_ar,
            dec1_1_ar, pligst_1_ar, strlig_1_ar, pheff_struc_ar)

        self.assert_all_values_in_array_within_range(
            tcflow_strucc1_ar, tcflow_strucc_1 - tolerance,
            tcflow_strucc_1 + tolerance, _IC_NODATA)

        insert_nodata_values_into_array(strmax_1_ar, _IC_NODATA)
        insert_nodata_values_into_array(dec1_1_ar, _IC_NODATA)
        insert_nodata_values_into_array(aminrl_1_ar, _SV_NODATA)

        tcflow_strucc1_ar = forage.calc_tcflow_strucc_1(
            aminrl_1_ar, aminrl_2_ar, strucc_1_ar, struce_1_1_ar,
            struce_1_2_ar, rnewas_1_1_ar, rnewas_2_1_ar, strmax_1_ar, defac_ar,
            dec1_1_ar, pligst_1_ar, strlig_1_ar, pheff_struc_ar)

        self.assert_all_values_in_array_within_range(
            tcflow_strucc1_ar, tcflow_strucc_1 - tolerance,
            tcflow_strucc_1 + tolerance, _IC_NODATA)

    def test_reclassify_nodata(self):
        """Test `reclassify_nodata`.

        Use the function `reclassify_nodata` to reset the nodata value of
        a raster.

        Raises:
            AssertionError if the nodata value of a raster following
                `reclassify_nodata` is not equal to the specified nodata type
            AssertionError if unique values in a raster contain more than the
                fill value and the correct nodata value

        Returns:
            None
        """
        from natcap.invest import forage

        fill_value = 0
        target_path = os.path.join(self.workspace_dir, 'target_raster.tif')
        create_random_raster(target_path, fill_value, fill_value)
        insert_nodata_values_into_raster(target_path, _TARGET_NODATA)

        new_nodata_value = -999
        forage.reclassify_nodata(target_path, new_nodata_value)
        result_nodata_value = pygeoprocessing.get_raster_info(
            target_path)['nodata'][0]
        self.assertEqual(
            new_nodata_value, result_nodata_value,
            msg="New nodata value does not match specified nodata value")

        # check unique values inside raster
        raster_values = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                target_path):
            raster_values.update(numpy.unique(raster_block))
        self.assertEqual(
            raster_values, set([fill_value, float(new_nodata_value)]),
            msg="Raster contains extraneous values")

        new_nodata_value = numpy.finfo('float32').min
        forage.reclassify_nodata(target_path, new_nodata_value)
        result_nodata_value = pygeoprocessing.get_raster_info(
            target_path)['nodata'][0]
        self.assertEqual(
            new_nodata_value, result_nodata_value,
            msg="New nodata value does not match specified nodata value")
        raster_values = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                target_path):
            raster_values.update(numpy.unique(raster_block))
        self.assertEqual(
            raster_values, set([fill_value, float(new_nodata_value)]),
            msg="Raster contains extraneous values")

        new_nodata_value = 8920
        forage.reclassify_nodata(target_path, new_nodata_value)
        insert_nodata_values_into_raster(target_path, new_nodata_value)
        result_nodata_value = pygeoprocessing.get_raster_info(
            target_path)['nodata'][0]
        self.assertEqual(
            new_nodata_value, result_nodata_value,
            msg="New nodata value does not match specified nodata value")
        raster_values = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                target_path):
            raster_values.update(numpy.unique(raster_block))
        self.assertEqual(
            raster_values, set([fill_value, float(new_nodata_value)]),
            msg="Raster contains extraneous values")

    def test_calc_respiration_mineral_flow(self):
        """Test `calc_respiration_mineral_flow`.

        Use the function `calc_respiration_mineral_flow` to calculate
        mineral flow of one element associated with respiration. Compare
        the result to values calculated by point-based verison defined
        here.

        Raises:
            AssertionError if `calc_respiration_mineral_flow` does not
                match values calculated by point-based version

        Returns:
            None
        """
        def respir_minr_flow_point(cflow, frac_co2, estatv, cstatv):
            co2_loss = cflow * frac_co2
            mineral_flow = co2_loss * estatv / cstatv
            return mineral_flow

        from natcap.invest import forage
        array_shape = (10, 10)
        tolerance = 0.0000000001

        # known values
        cflow = 15.2006601
        frac_co2 = 0.0146
        estatv = 0.7776
        cstatv = 155.5253
        mineral_flow = respir_minr_flow_point(cflow, frac_co2, estatv, cstatv)

        cflow_ar = numpy.full(array_shape, cflow)
        frac_co2_ar = numpy.full(array_shape, frac_co2)
        estatv_ar = numpy.full(array_shape, estatv)
        cstatv_ar = numpy.full(array_shape, cstatv)
        mineral_flow_ar = forage.calc_respiration_mineral_flow(
            cflow_ar, frac_co2_ar, estatv_ar, cstatv_ar)

        self.assert_all_values_in_array_within_range(
            mineral_flow_ar, mineral_flow - tolerance,
            mineral_flow + tolerance, _IC_NODATA)

        insert_nodata_values_into_array(cflow_ar, _IC_NODATA)
        insert_nodata_values_into_array(frac_co2_ar, _IC_NODATA)
        insert_nodata_values_into_array(estatv_ar, _SV_NODATA)
        insert_nodata_values_into_array(cstatv_ar, _SV_NODATA)

        mineral_flow_ar = forage.calc_respiration_mineral_flow(
            cflow_ar, frac_co2_ar, estatv_ar, cstatv_ar)

        self.assert_all_values_in_array_within_range(
            mineral_flow_ar, mineral_flow - tolerance,
            mineral_flow + tolerance, _IC_NODATA)

    def test_update_gross_mineralization(self):
        """Test `update_gross_mineralization`.

        Test the function `update_gross_mineralization`  against values
        calculated by hand.

        Raises:
            AssertionError if `update_gross_mineralization` does not match
                values calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        array_shape = (10, 10)
        tolerance = 0.0000001

        # known values
        gross_mineralization = 0.0481
        mineral_flow = 0.00044
        gromin_updated = gross_mineralization + mineral_flow

        gross_mineralization_ar = numpy.full(array_shape, gross_mineralization)
        mineral_flow_ar = numpy.full(array_shape, mineral_flow)

        gromin_updated_ar = forage.update_gross_mineralization(
            gross_mineralization_ar, mineral_flow_ar)
        self.assert_all_values_in_array_within_range(
            gromin_updated_ar, gromin_updated - tolerance,
            gromin_updated + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(
            gross_mineralization_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(mineral_flow_ar, _IC_NODATA)

        gromin_updated_ar = forage.update_gross_mineralization(
            gross_mineralization_ar, mineral_flow_ar)
        self.assert_all_values_in_array_within_range(
            gromin_updated_ar, gromin_updated - tolerance,
            gromin_updated + tolerance, _TARGET_NODATA)

        # known values
        gross_mineralization = 0.048881
        mineral_flow = -0.00674
        gromin_updated = gross_mineralization

        gross_mineralization_ar = numpy.full(array_shape, gross_mineralization)
        mineral_flow_ar = numpy.full(array_shape, mineral_flow)

        gromin_updated_ar = forage.update_gross_mineralization(
            gross_mineralization_ar, mineral_flow_ar)
        self.assert_all_values_in_array_within_range(
            gromin_updated_ar, gromin_updated - tolerance,
            gromin_updated + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(
            gross_mineralization_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(mineral_flow_ar, _IC_NODATA)

        gromin_updated_ar = forage.update_gross_mineralization(
            gross_mineralization_ar, mineral_flow_ar)
        self.assert_all_values_in_array_within_range(
            gromin_updated_ar, gromin_updated - tolerance,
            gromin_updated + tolerance, _TARGET_NODATA)

    def test_calc_net_cflow(self):
        """Test `calc_net_cflow`.

        Test `calc_net_cflow` against value calculated by hand.

        Raises:
            AssertionError if `calc_net_cflow` does not match values
                calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        array_shape = (10, 10)
        tolerance = 0.0000001

        # known values
        cflow = 5.2006601
        frac_co2 = 0.0182
        net_cflow = cflow - (cflow * frac_co2)

        cflow_ar = numpy.full(array_shape, cflow)
        frac_co2_ar = numpy.full(array_shape, frac_co2)
        net_cflow_ar = forage.calc_net_cflow(cflow_ar, frac_co2_ar)

        self.assert_all_values_in_array_within_range(
            net_cflow_ar, net_cflow - tolerance, net_cflow + tolerance,
            _IC_NODATA)

        insert_nodata_values_into_array(cflow_ar, _IC_NODATA)
        insert_nodata_values_into_array(frac_co2_ar, _IC_NODATA)

        net_cflow_ar = forage.calc_net_cflow(cflow_ar, frac_co2_ar)

        self.assert_all_values_in_array_within_range(
            net_cflow_ar, net_cflow - tolerance, net_cflow + tolerance,
            _IC_NODATA)

    def test_calc_tcflow_surface(self):
        """Test `calc_tcflow_surface`.

        Test `calc_tcflow_surface` against value calculated by point-based
        version.

        Raises:
            AssertionError if `calc_tcflow_surface` does not match value
                calculated by point-based version

        Returns:
            None
        """
        def calc_tcflow_surface_point(
                aminrl_1, aminrl_2, metabc_1, metabe_1_1, metabe_1_2,
                rceto1_1, rceto1_2, defac, dec2_1, pheff_metab):
            """Point implementation of `calc_tcflow_surface`."""
            decompose_mask = (
                ((aminrl_1 > 0.0000001) | (
                    (metabc_1 / metabe_1_1) <= rceto1_1)) &
                ((aminrl_2 > 0.0000001) | (
                    (metabc_1 / metabe_1_2) <= rceto1_2)))  # line 194 Litdec.f
            if decompose_mask:
                tcflow_metabc_1 = numpy.clip(
                    (metabc_1 * defac * dec2_1 * 0.020833 * pheff_metab), 0,
                    metabc_1)
            else:
                tcflow_metabc_1 = 0.
            return tcflow_metabc_1
        from natcap.invest import forage
        array_shape = (10, 10)
        tolerance = 0.00001

        # known values, decomposition can occur
        aminrl_1 = 5.8821
        aminrl_2 = 0.04781
        metabc_1 = 169.22
        metabe_1_1 = 0.7776
        metabe_1_2 = 0.3111
        rceto1_1 = 5.29
        rceto1_2 = 2.92
        defac = 0.822
        dec2_1 = 3.9
        pheff_metab = 0.9917

        tcflow_metabc_1_point = calc_tcflow_surface_point(
            aminrl_1, aminrl_2, metabc_1, metabe_1_1, metabe_1_2,
            rceto1_1, rceto1_2, defac, dec2_1, pheff_metab)

        # raster inputs
        aminrl_1_ar = numpy.full(array_shape, aminrl_1)
        aminrl_2_ar = numpy.full(array_shape, aminrl_2)
        metabc_1_ar = numpy.full(array_shape, metabc_1)
        metabe_1_1_ar = numpy.full(array_shape, metabe_1_1)
        metabe_1_2_ar = numpy.full(array_shape, metabe_1_2)
        rceto1_1_ar = numpy.full(array_shape, rceto1_1)
        rceto1_2_ar = numpy.full(array_shape, rceto1_2)
        defac_ar = numpy.full(array_shape, defac)
        dec2_1_ar = numpy.full(array_shape, dec2_1)
        pheff_metab_ar = numpy.full(array_shape, pheff_metab)

        tcflow_metabc_1_ar = forage.calc_tcflow_surface(
            aminrl_1_ar, aminrl_2_ar, metabc_1_ar, metabe_1_1_ar,
            metabe_1_2_ar, rceto1_1_ar, rceto1_2_ar, defac_ar, dec2_1_ar,
            pheff_metab_ar)

        self.assert_all_values_in_array_within_range(
            tcflow_metabc_1_ar, tcflow_metabc_1_point - tolerance,
            tcflow_metabc_1_point + tolerance, _IC_NODATA)

        insert_nodata_values_into_array(aminrl_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(defac_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(metabe_1_2_ar, _SV_NODATA)
        insert_nodata_values_into_array(metabe_1_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(metabc_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(pheff_metab_ar, _TARGET_NODATA)

        tcflow_metabc_1_ar = forage.calc_tcflow_surface(
            aminrl_1_ar, aminrl_2_ar, metabc_1_ar, metabe_1_1_ar,
            metabe_1_2_ar, rceto1_1_ar, rceto1_2_ar, defac_ar, dec2_1_ar,
            pheff_metab_ar)

        self.assert_all_values_in_array_within_range(
            tcflow_metabc_1_ar, tcflow_metabc_1_point - tolerance,
            tcflow_metabc_1_point + tolerance, _IC_NODATA)

        # known values, no decomposition
        aminrl_1 = 0.
        aminrl_2 = 0.
        metabc_1 = 169.22
        metabe_1_1 = 0.7776
        metabe_1_2 = 0.3111
        rceto1_1 = 200.
        rceto1_2 = 400.
        defac = 0.822
        dec2_1 = 3.9
        pheff_metab = 0.9917

        tcflow_metabc_1_point = calc_tcflow_surface_point(
            aminrl_1, aminrl_2, metabc_1, metabe_1_1, metabe_1_2,
            rceto1_1, rceto1_2, defac, dec2_1, pheff_metab)

        # raster inputs
        aminrl_1_ar = numpy.full(array_shape, aminrl_1)
        aminrl_2_ar = numpy.full(array_shape, aminrl_2)
        metabc_1_ar = numpy.full(array_shape, metabc_1)
        metabe_1_1_ar = numpy.full(array_shape, metabe_1_1)
        metabe_1_2_ar = numpy.full(array_shape, metabe_1_2)
        rceto1_1_ar = numpy.full(array_shape, rceto1_1)
        rceto1_2_ar = numpy.full(array_shape, rceto1_2)
        defac_ar = numpy.full(array_shape, defac)
        dec2_1_ar = numpy.full(array_shape, dec2_1)
        pheff_metab_ar = numpy.full(array_shape, pheff_metab)

        tcflow_metabc_1_ar = forage.calc_tcflow_surface(
            aminrl_1_ar, aminrl_2_ar, metabc_1_ar, metabe_1_1_ar,
            metabe_1_2_ar, rceto1_1_ar, rceto1_2_ar, defac_ar, dec2_1_ar,
            pheff_metab_ar)

        self.assert_all_values_in_array_within_range(
            tcflow_metabc_1_ar, tcflow_metabc_1_point - tolerance,
            tcflow_metabc_1_point + tolerance, _IC_NODATA)

        insert_nodata_values_into_array(aminrl_2_ar, _SV_NODATA)
        insert_nodata_values_into_array(defac_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(rceto1_2_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(metabe_1_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(dec2_1_ar, _IC_NODATA)
        insert_nodata_values_into_array(pheff_metab_ar, _TARGET_NODATA)

        tcflow_metabc_1_ar = forage.calc_tcflow_surface(
            aminrl_1_ar, aminrl_2_ar, metabc_1_ar, metabe_1_1_ar,
            metabe_1_2_ar, rceto1_1_ar, rceto1_2_ar, defac_ar, dec2_1_ar,
            pheff_metab_ar)

        self.assert_all_values_in_array_within_range(
            tcflow_metabc_1_ar, tcflow_metabc_1_point - tolerance,
            tcflow_metabc_1_point + tolerance, _IC_NODATA)

    def test_calc_tcflow_soil(self):
        """Test `calc_tcflow_soil`.

        Test `calc_tcflow_soil` against value calculated by point-based
        version.

        Raises:
            AssertionError if `calc_tcflow_soil` does not match value
                calculated by point-based version

        Returns:
            None
        """
        def calc_tcflow_soil_point(
                aminrl_1, aminrl_2, metabc_2, metabe_2_1, metabe_2_2, rceto1_1,
                rceto1_2, defac, dec2_2, pheff_metab, anerb):
            """Point implementation of `calc_tcflow_soil`."""
            decompose_mask = (
                ((aminrl_1 > 0.0000001) | (
                    (metabc_2 / metabe_2_1) <= rceto1_1)) &
                ((aminrl_2 > 0.0000001) | (
                    (metabc_2 / metabe_2_2) <= rceto1_2)))  # line 194 Litdec.f
            if decompose_mask:
                tcflow_metabc_2 = numpy.clip(
                    (metabc_2 * defac * dec2_2 * 0.020833 * pheff_metab *
                        anerb), 0, metabc_2)
            else:
                tcflow_metabc_2 = 0.
            return tcflow_metabc_2
        from natcap.invest import forage
        array_shape = (10, 10)
        tolerance = 0.00001

        # known values, decomposition can occur
        aminrl_1 = 5.8821
        aminrl_2 = 0.04781
        metabc_2 = 169.22
        metabe_2_1 = 0.7776
        metabe_2_2 = 0.3111
        rceto1_1 = 5.29
        rceto1_2 = 2.92
        defac = 0.822
        dec2_2 = 3.9
        pheff_metab = 0.9917
        anerb = 0.3

        tcflow_metabc_2_point = calc_tcflow_soil_point(
            aminrl_1, aminrl_2, metabc_2, metabe_2_1, metabe_2_2,
            rceto1_1, rceto1_2, defac, dec2_2, pheff_metab, anerb)

        # raster inputs
        aminrl_1_ar = numpy.full(array_shape, aminrl_1)
        aminrl_2_ar = numpy.full(array_shape, aminrl_2)
        metabc_2_ar = numpy.full(array_shape, metabc_2)
        metabe_2_1_ar = numpy.full(array_shape, metabe_2_1)
        metabe_2_2_ar = numpy.full(array_shape, metabe_2_2)
        rceto1_1_ar = numpy.full(array_shape, rceto1_1)
        rceto1_2_ar = numpy.full(array_shape, rceto1_2)
        defac_ar = numpy.full(array_shape, defac)
        dec2_2_ar = numpy.full(array_shape, dec2_2)
        pheff_metab_ar = numpy.full(array_shape, pheff_metab)
        anerb_ar = numpy.full(array_shape, anerb)

        tcflow_metabc_2_ar = forage.calc_tcflow_soil(
            aminrl_1_ar, aminrl_2_ar, metabc_2_ar, metabe_2_1_ar,
            metabe_2_2_ar, rceto1_1_ar, rceto1_2_ar, defac_ar, dec2_2_ar,
            pheff_metab_ar, anerb_ar)
        self.assert_all_values_in_array_within_range(
            tcflow_metabc_2_ar, tcflow_metabc_2_point - tolerance,
            tcflow_metabc_2_point + tolerance, _IC_NODATA)

        insert_nodata_values_into_array(aminrl_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(defac_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(metabe_2_2_ar, _SV_NODATA)
        insert_nodata_values_into_array(metabe_2_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(anerb_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(metabc_2_ar, _SV_NODATA)
        insert_nodata_values_into_array(pheff_metab_ar, _TARGET_NODATA)

        tcflow_metabc_2_ar = forage.calc_tcflow_soil(
            aminrl_1_ar, aminrl_2_ar, metabc_2_ar, metabe_2_1_ar,
            metabe_2_2_ar, rceto1_1_ar, rceto1_2_ar, defac_ar, dec2_2_ar,
            pheff_metab_ar, anerb_ar)
        self.assert_all_values_in_array_within_range(
            tcflow_metabc_2_ar, tcflow_metabc_2_point - tolerance,
            tcflow_metabc_2_point + tolerance, _IC_NODATA)

        # known values, no decomposition
        aminrl_1 = 0.
        aminrl_2 = 0.
        metabc_2 = 169.22
        metabe_2_1 = 0.7776
        metabe_2_2 = 0.3111
        rceto1_1 = 200.
        rceto1_2 = 400.
        defac = 0.822
        dec2_2 = 3.9
        pheff_metab = 0.9917

        tcflow_metabc_2_point = calc_tcflow_soil_point(
            aminrl_1, aminrl_2, metabc_2, metabe_2_1, metabe_2_2,
            rceto1_1, rceto1_2, defac, dec2_2, pheff_metab, anerb)

        # raster inputs
        aminrl_1_ar = numpy.full(array_shape, aminrl_1)
        aminrl_2_ar = numpy.full(array_shape, aminrl_2)
        metabc_2_ar = numpy.full(array_shape, metabc_2)
        metabe_2_1_ar = numpy.full(array_shape, metabe_2_1)
        metabe_2_2_ar = numpy.full(array_shape, metabe_2_2)
        rceto1_1_ar = numpy.full(array_shape, rceto1_1)
        rceto1_2_ar = numpy.full(array_shape, rceto1_2)
        defac_ar = numpy.full(array_shape, defac)
        dec2_2_ar = numpy.full(array_shape, dec2_2)
        pheff_metab_ar = numpy.full(array_shape, pheff_metab)
        anerb_ar = numpy.full(array_shape, anerb)

        tcflow_metabc_2_ar = forage.calc_tcflow_soil(
            aminrl_1_ar, aminrl_2_ar, metabc_2_ar, metabe_2_1_ar,
            metabe_2_2_ar, rceto1_1_ar, rceto1_2_ar, defac_ar, dec2_2_ar,
            pheff_metab_ar, anerb_ar)
        self.assert_all_values_in_array_within_range(
            tcflow_metabc_2_ar, tcflow_metabc_2_point - tolerance,
            tcflow_metabc_2_point + tolerance, _IC_NODATA)

        insert_nodata_values_into_array(aminrl_2_ar, _SV_NODATA)
        insert_nodata_values_into_array(defac_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(rceto1_2_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(metabe_2_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(dec2_2_ar, _IC_NODATA)
        insert_nodata_values_into_array(pheff_metab_ar, _TARGET_NODATA)

        tcflow_metabc_2_ar = forage.calc_tcflow_soil(
            aminrl_1_ar, aminrl_2_ar, metabc_2_ar, metabe_2_1_ar,
            metabe_2_2_ar, rceto1_1_ar, rceto1_2_ar, defac_ar, dec2_2_ar,
            pheff_metab_ar, anerb_ar)
        self.assert_all_values_in_array_within_range(
            tcflow_metabc_2_ar, tcflow_metabc_2_point - tolerance,
            tcflow_metabc_2_point + tolerance, _IC_NODATA)

    def test_belowground_ratio(self):
        """Test `_belowground_ratio`.

        Use the function `_belowground_ratio` to calculate required C/iel
        ratio of belowground decomposing material.  Compare calculate values
        to values calculated by point-based version.

        Raises:
            AssertionError if `_belowground_ratio` does not match values
                calculated by point-based version

        Returns:
            None
        """
        from natcap.invest import forage
        array_shape = (10, 10)
        tolerance = 0.00001

        # known values, aminrl > varat_3_iel
        aminrl = 5.928
        varat_1_iel = 14.
        varat_2_iel = 3.
        varat_3_iel = 2.

        belowground_point = bgdrat_point(
            aminrl, varat_1_iel, varat_2_iel, varat_3_iel)

        # array inputs
        aminrl_ar = numpy.full(array_shape, aminrl)
        varat_1_iel_ar = numpy.full(array_shape, varat_1_iel)
        varat_2_iel_ar = numpy.full(array_shape, varat_2_iel)
        varat_3_iel_ar = numpy.full(array_shape, varat_3_iel)

        belowground_ratio = forage._belowground_ratio(
            aminrl_ar, varat_1_iel_ar, varat_2_iel_ar, varat_3_iel_ar)
        self.assert_all_values_in_array_within_range(
            belowground_ratio, belowground_point - tolerance,
            belowground_point + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(aminrl_ar, _SV_NODATA)
        insert_nodata_values_into_array(varat_1_iel_ar, _IC_NODATA)
        insert_nodata_values_into_array(varat_2_iel_ar, _IC_NODATA)
        insert_nodata_values_into_array(varat_3_iel_ar, _IC_NODATA)

        belowground_ratio = forage._belowground_ratio(
            aminrl_ar, varat_1_iel_ar, varat_2_iel_ar, varat_3_iel_ar)
        self.assert_all_values_in_array_within_range(
            belowground_ratio, belowground_point - tolerance,
            belowground_point + tolerance, _TARGET_NODATA)

        # no mineral source
        aminrl = 0.
        varat_1_iel = 14.
        varat_2_iel = 3.
        varat_3_iel = 2.

        belowground_point = bgdrat_point(
            aminrl, varat_1_iel, varat_2_iel, varat_3_iel)

        # array inputs
        aminrl_ar = numpy.full(array_shape, aminrl)
        varat_1_iel_ar = numpy.full(array_shape, varat_1_iel)
        varat_2_iel_ar = numpy.full(array_shape, varat_2_iel)
        varat_3_iel_ar = numpy.full(array_shape, varat_3_iel)

        belowground_ratio = forage._belowground_ratio(
            aminrl_ar, varat_1_iel_ar, varat_2_iel_ar, varat_3_iel_ar)
        self.assert_all_values_in_array_within_range(
            belowground_ratio, belowground_point - tolerance,
            belowground_point + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(aminrl_ar, _SV_NODATA)
        insert_nodata_values_into_array(varat_1_iel_ar, _IC_NODATA)
        insert_nodata_values_into_array(varat_2_iel_ar, _IC_NODATA)
        insert_nodata_values_into_array(varat_3_iel_ar, _IC_NODATA)

        belowground_ratio = forage._belowground_ratio(
            aminrl_ar, varat_1_iel_ar, varat_2_iel_ar, varat_3_iel_ar)
        self.assert_all_values_in_array_within_range(
            belowground_ratio, belowground_point - tolerance,
            belowground_point + tolerance, _TARGET_NODATA)

        # known values, aminrl < varat_3_iel
        aminrl = 1.9917
        varat_1_iel = 14.
        varat_2_iel = 5.
        varat_3_iel = 3.

        belowground_point = bgdrat_point(
            aminrl, varat_1_iel, varat_2_iel, varat_3_iel)

        # array inputs
        aminrl_ar = numpy.full(array_shape, aminrl)
        varat_1_iel_ar = numpy.full(array_shape, varat_1_iel)
        varat_2_iel_ar = numpy.full(array_shape, varat_2_iel)
        varat_3_iel_ar = numpy.full(array_shape, varat_3_iel)

        belowground_ratio = forage._belowground_ratio(
            aminrl_ar, varat_1_iel_ar, varat_2_iel_ar, varat_3_iel_ar)
        self.assert_all_values_in_array_within_range(
            belowground_ratio, belowground_point - tolerance,
            belowground_point + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(aminrl_ar, _SV_NODATA)
        insert_nodata_values_into_array(varat_1_iel_ar, _IC_NODATA)
        insert_nodata_values_into_array(varat_2_iel_ar, _IC_NODATA)
        insert_nodata_values_into_array(varat_3_iel_ar, _IC_NODATA)

        belowground_ratio = forage._belowground_ratio(
            aminrl_ar, varat_1_iel_ar, varat_2_iel_ar, varat_3_iel_ar)
        self.assert_all_values_in_array_within_range(
            belowground_ratio, belowground_point - tolerance,
            belowground_point + tolerance, _TARGET_NODATA)

    def test_calc_surface_som2_ratio(self):
        """Test `calc_surface_som2_ratio`.

        Use the function `calc_surface_som2_ratio` to calculate the required
        ratio for material entering surface SOM2. Compare the calculated
        value to value calculated by hand.

        Raises:
            AssertionError if `calc_surface_som2_ratio` does not match value
                calculated by hand

        Returns:
            None
        """
        from natcap.invest import forage
        array_shape = (10, 10)
        tolerance = 0.00001

        # known values, calc term > rad1p_3
        som1c_1 = 12.8192
        som1e_1 = 1.5752
        rad1p_1 = 12.
        rad1p_2 = 3.
        rad1p_3 = 5.
        pcemic1_2 = 10.
        rceto2_surface = 14.552565

        # array inputs
        som1c_1_ar = numpy.full(array_shape, som1c_1)
        som1e_1_ar = numpy.full(array_shape, som1e_1)
        rad1p_1_ar = numpy.full(array_shape, rad1p_1)
        rad1p_2_ar = numpy.full(array_shape, rad1p_2)
        rad1p_3_ar = numpy.full(array_shape, rad1p_3)
        pcemic1_2_ar = numpy.full(array_shape, pcemic1_2)

        receto2_surface_ar = forage.calc_surface_som2_ratio(
            som1c_1_ar,  som1e_1_ar, rad1p_1_ar, rad1p_2_ar, rad1p_3_ar,
            pcemic1_2_ar)
        self.assert_all_values_in_array_within_range(
            receto2_surface_ar, rceto2_surface - tolerance,
            rceto2_surface + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(som1c_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(som1e_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(rad1p_1_ar, _IC_NODATA)
        insert_nodata_values_into_array(rad1p_2_ar, _IC_NODATA)
        insert_nodata_values_into_array(rad1p_3_ar, _IC_NODATA)
        insert_nodata_values_into_array(pcemic1_2_ar, _IC_NODATA)

        receto2_surface_ar = forage.calc_surface_som2_ratio(
            som1c_1_ar,  som1e_1_ar, rad1p_1_ar, rad1p_2_ar, rad1p_3_ar,
            pcemic1_2_ar)
        self.assert_all_values_in_array_within_range(
            receto2_surface_ar, rceto2_surface - tolerance,
            rceto2_surface + tolerance, _TARGET_NODATA)

        # known values, calc term < rad1p_3
        som1c_1 = 12.8192
        som1e_1 = 2.
        rad1p_1 = 8.
        rad1p_2 = 3.
        rad1p_3 = 5.
        pcemic1_2 = 10.
        rceto2_surface = 5.

        # array inputs
        som1c_1_ar = numpy.full(array_shape, som1c_1)
        som1e_1_ar = numpy.full(array_shape, som1e_1)
        rad1p_1_ar = numpy.full(array_shape, rad1p_1)
        rad1p_2_ar = numpy.full(array_shape, rad1p_2)
        rad1p_3_ar = numpy.full(array_shape, rad1p_3)
        pcemic1_2_ar = numpy.full(array_shape, pcemic1_2)

        receto2_surface_ar = forage.calc_surface_som2_ratio(
            som1c_1_ar,  som1e_1_ar, rad1p_1_ar, rad1p_2_ar, rad1p_3_ar,
            pcemic1_2_ar)
        self.assert_all_values_in_array_within_range(
            receto2_surface_ar, rceto2_surface - tolerance,
            rceto2_surface + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(som1c_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(som1e_1_ar, _SV_NODATA)
        insert_nodata_values_into_array(rad1p_1_ar, _IC_NODATA)
        insert_nodata_values_into_array(rad1p_2_ar, _IC_NODATA)
        insert_nodata_values_into_array(rad1p_3_ar, _IC_NODATA)
        insert_nodata_values_into_array(pcemic1_2_ar, _IC_NODATA)

        receto2_surface_ar = forage.calc_surface_som2_ratio(
            som1c_1_ar,  som1e_1_ar, rad1p_1_ar, rad1p_2_ar, rad1p_3_ar,
            pcemic1_2_ar)
        self.assert_all_values_in_array_within_range(
            receto2_surface_ar, rceto2_surface - tolerance,
            rceto2_surface + tolerance, _TARGET_NODATA)

    def test_calc_c_leach(self):
        """Test `calc_c_leach`.

        Use the function `calc_c_leach` to calculate the C leaching
        from soil SOM1 into stream flow during decomposition.  Compare
        calculated values to value calculated by hand.

        Raises:
            AssertionError if `calc_c_leach` does not match value
                calculated by hand.

        Returns:
            None
        """
        from natcap.invest import forage
        array_shape = (10, 10)
        tolerance = 0.00001

        # known values, linten > 1
        amov_2 = 63.1
        tcflow = 40.38
        omlech_3 = 60.
        orglch = 0.07
        cleach = 2.8266

        # array inputs
        amov_2_ar = numpy.full(array_shape, amov_2)
        tcflow_ar = numpy.full(array_shape, tcflow)
        omlech_3_ar = numpy.full(array_shape, omlech_3)
        orglch_ar = numpy.full(array_shape, orglch)

        cleach_ar = forage.calc_c_leach(
            amov_2_ar, tcflow_ar, omlech_3_ar, orglch_ar)
        self.assert_all_values_in_array_within_range(
            cleach_ar, cleach - tolerance, cleach + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(amov_2_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(tcflow_ar, _IC_NODATA)
        insert_nodata_values_into_array(omlech_3_ar, _IC_NODATA)
        insert_nodata_values_into_array(orglch_ar, _IC_NODATA)

        cleach_ar = forage.calc_c_leach(
            amov_2_ar, tcflow_ar, omlech_3_ar, orglch_ar)
        self.assert_all_values_in_array_within_range(
            cleach_ar, cleach - tolerance, cleach + tolerance, _TARGET_NODATA)

        # known values, linten < 1
        amov_2 = 10.5
        tcflow = 40.38
        omlech_3 = 60.
        orglch = 0.07
        cleach = 0.494655

        # array inputs
        amov_2_ar = numpy.full(array_shape, amov_2)
        tcflow_ar = numpy.full(array_shape, tcflow)
        omlech_3_ar = numpy.full(array_shape, omlech_3)
        orglch_ar = numpy.full(array_shape, orglch)

        cleach_ar = forage.calc_c_leach(
            amov_2_ar, tcflow_ar, omlech_3_ar, orglch_ar)
        self.assert_all_values_in_array_within_range(
            cleach_ar, cleach - tolerance, cleach + tolerance, _TARGET_NODATA)

        insert_nodata_values_into_array(amov_2_ar, _TARGET_NODATA)
        insert_nodata_values_into_array(tcflow_ar, _IC_NODATA)
        insert_nodata_values_into_array(omlech_3_ar, _IC_NODATA)
        insert_nodata_values_into_array(orglch_ar, _IC_NODATA)

        cleach_ar = forage.calc_c_leach(
            amov_2_ar, tcflow_ar, omlech_3_ar, orglch_ar)
        self.assert_all_values_in_array_within_range(
            cleach_ar, cleach - tolerance, cleach + tolerance, _TARGET_NODATA)

    def test_remove_leached_iel(self):
        """Test `remove_leached_iel`.

        Use the function `remove_leached_iel` to remove leached nutrients from
        SOM1 during decomposition. Test the calculated values against values
        calculated by hand.

        Raises:
            AssertionError if calculated values do not match values calculated
                by hand
        """
        from natcap.invest import forage
        nrows = 10
        ncols = 10
        tolerance = 0.00001

        # known values, leaching N
        som1c_2 = 10.5
        som1e_2_iel = 40.38
        cleach = 0.494655
        iel = 1
        d_som1e_2_iel_before = 50.22
        d_som1e_2_iel_after = 49.2688491

        # raster inputs
        som1c_2_path = os.path.join(self.workspace_dir, 'som1c_2.tif')
        som1e_2_iel_path = os.path.join(self.workspace_dir, 'som1e_2_iel.tif')
        cleach_path = os.path.join(self.workspace_dir, 'cleach.tif')
        d_som1e_2_iel_path = os.path.join(
            self.workspace_dir, 'd_som1e_2_iel.tif')
        create_random_raster(som1c_2_path, som1c_2, som1c_2)
        create_random_raster(som1e_2_iel_path, som1e_2_iel, som1e_2_iel)
        create_random_raster(cleach_path, cleach, cleach)
        create_random_raster(
            d_som1e_2_iel_path, d_som1e_2_iel_before, d_som1e_2_iel_before)

        forage.remove_leached_iel(
            som1c_2_path, som1e_2_iel_path, cleach_path, d_som1e_2_iel_path,
            iel)
        self.assert_all_values_in_raster_within_range(
            d_som1e_2_iel_path, d_som1e_2_iel_after - tolerance,
            d_som1e_2_iel_after + tolerance, _IC_NODATA)

        create_random_raster(
            d_som1e_2_iel_path, d_som1e_2_iel_before, d_som1e_2_iel_before)
        insert_nodata_values_into_raster(som1c_2_path, _SV_NODATA)
        insert_nodata_values_into_raster(som1e_2_iel_path, _SV_NODATA)
        insert_nodata_values_into_raster(cleach_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(d_som1e_2_iel_path, _IC_NODATA)

        forage.remove_leached_iel(
            som1c_2_path, som1e_2_iel_path, cleach_path, d_som1e_2_iel_path,
            iel)
        self.assert_all_values_in_raster_within_range(
            d_som1e_2_iel_path, d_som1e_2_iel_after - tolerance,
            d_som1e_2_iel_after + tolerance, _IC_NODATA)

        # known values, leaching P
        som1c_2 = 10.5
        som1e_2_iel = 40.38
        cleach = 0.494655
        iel = 2
        d_som1e_2_iel_before = 50.22
        d_som1e_2_iel_after = 50.16564852

        # raster inputs
        som1c_2_path = os.path.join(self.workspace_dir, 'som1c_2.tif')
        som1e_2_iel_path = os.path.join(self.workspace_dir, 'som1e_2_iel.tif')
        cleach_path = os.path.join(self.workspace_dir, 'cleach.tif')
        d_som1e_2_iel_path = os.path.join(
            self.workspace_dir, 'd_som1e_2_iel.tif')
        create_random_raster(som1c_2_path, som1c_2, som1c_2)
        create_random_raster(som1e_2_iel_path, som1e_2_iel, som1e_2_iel)
        create_random_raster(cleach_path, cleach, cleach)
        create_random_raster(
            d_som1e_2_iel_path, d_som1e_2_iel_before, d_som1e_2_iel_before)

        forage.remove_leached_iel(
            som1c_2_path, som1e_2_iel_path, cleach_path, d_som1e_2_iel_path,
            iel)
        self.assert_all_values_in_raster_within_range(
            d_som1e_2_iel_path, d_som1e_2_iel_after - tolerance,
            d_som1e_2_iel_after + tolerance, _IC_NODATA)

        create_random_raster(
            d_som1e_2_iel_path, d_som1e_2_iel_before, d_som1e_2_iel_before)
        insert_nodata_values_into_raster(som1c_2_path, _SV_NODATA)
        insert_nodata_values_into_raster(som1e_2_iel_path, _SV_NODATA)
        insert_nodata_values_into_raster(cleach_path, _TARGET_NODATA)
        insert_nodata_values_into_raster(d_som1e_2_iel_path, _IC_NODATA)

        forage.remove_leached_iel(
            som1c_2_path, som1e_2_iel_path, cleach_path, d_som1e_2_iel_path,
            iel)
        self.assert_all_values_in_raster_within_range(
            d_som1e_2_iel_path, d_som1e_2_iel_after - tolerance,
            d_som1e_2_iel_after + tolerance, _IC_NODATA)
