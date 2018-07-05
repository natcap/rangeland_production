"""InVEST forage model tests."""

import unittest
import tempfile
import shutil
import os
import sys

import numpy
from osgeo import osr
from osgeo import gdal

import pygeoprocessing

SAMPLE_DATA = "C:/Users/ginge/Documents/NatCap/sample_inputs"
REGRESSION_DATA = "C:/Users/ginge/Documents/NatCap/regression_test_data"

_TARGET_NODATA = -1.0
_IC_NODATA = numpy.finfo('float32').min


class foragetests(unittest.TestCase):
    """Regression tests for InVEST forage model."""

    def setUp(self):
        """Create temporary workspace directory."""
        self.workspace_dir = tempfile.mkdtemp()

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

    def test_model_runs(self):
        """Launch forage model, ensure it runs!."""
        from natcap.invest import forage

        if not os.path.exists(SAMPLE_DATA):
            self.fail(
                "Sample input directory not found at %s" % SAMPLE_DATA)

        args = foragetests.generate_base_args(self.workspace_dir)
        # forage.execute(args)

    def test_shortwave_radiation(self):
        """Test calculation of shortwave radiation outside the atmosphere.

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

        # create_test_raster of 1 pixel
        geotransform = [0, 1, 0, 44.5, 0, 1]
        n_cols = 1
        n_rows = 1
        n_bands = 1
        datatype = gdal.GDT_Float32
        projection = osr.SpatialReference()
        projection.SetWellKnownGeogCS('WGS84')
        driver = gdal.GetDriverByName('GTiff')
        target_raster = driver.Create(
            template_raster.encode('utf-8'), n_cols, n_rows, n_bands, datatype)
        target_raster.SetProjection(projection.ExportToWkt())
        target_raster.SetGeoTransform(geotransform)
        target_band = target_raster.GetRasterBand(1)
        target_band.SetNoDataValue(_TARGET_NODATA)
        target_band.Fill(fill_value)
        target_raster = None

        month = 5
        shwave_path = os.path.join(self.workspace_dir, 'shwave.tif')
        forage._shortwave_radiation(template_raster, month, shwave_path)

        # assert the value in the raster `shwave_path` is equal to value
        # calculated by hand
        result_set = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                shwave_path):
            result_set.update(numpy.unique(raster_block))
        assert len(result_set) == 1, (
            "One unique value expected in shortwave radiation raster")
        test_result = list(result_set)[0]
        assert abs(test_result - 990.7401) < 0.01, (
            "Test result does not match expected value")

    def test_calc_ompc(self):
        """Test the function _calc_ompc against value calculated by hand.

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

        def create_test_raster(target_path, fill_value):
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

        som1c_2_path = os.path.join(self.workspace_dir, 'som1c_2.tif')
        som2c_2_path = os.path.join(self.workspace_dir, 'som2c_2.tif')
        som3c_path = os.path.join(self.workspace_dir, 'som3c.tif')
        bulk_d_path = os.path.join(self.workspace_dir, 'bulkd.tif')
        edepth_path = os.path.join(self.workspace_dir, 'edepth.tif')

        create_test_raster(som1c_2_path, 42.109)
        create_test_raster(som2c_2_path, 959.1091)
        create_test_raster(som3c_path, 588.0574)
        create_test_raster(bulk_d_path, 1.5)
        create_test_raster(edepth_path, 0.2)

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
        assert len(result_set) == 1, (
            "One unique value expected in organic matter raster")
        test_result = list(result_set)[0]
        assert abs(test_result - 0.913304) < 0.0001, (
            "Test result does not match expected value")

    def test_calc_afiel(self):
        """Test the function _calc_afiel against value calculated by hand.

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

        def create_test_raster(target_path, fill_value):
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

        sand_path = os.path.join(self.workspace_dir, 'sand.tif')
        silt_path = os.path.join(self.workspace_dir, 'silt.tif')
        clay_path = os.path.join(self.workspace_dir, 'clay.tif')
        ompc_path = os.path.join(self.workspace_dir, 'ompc.tif')
        bulkd_path = os.path.join(self.workspace_dir, 'bulkd.tif')

        create_test_raster(sand_path, 0.39)
        create_test_raster(silt_path, 0.41)
        create_test_raster(clay_path, 0.2)
        create_test_raster(ompc_path, 0.913304)
        create_test_raster(bulkd_path, 1.5)

        afiel_path = os.path.join(self.workspace_dir, 'afiel.tif')

        forage._calc_afiel(
            sand_path, silt_path, clay_path, ompc_path, bulkd_path, afiel_path)

        # assert the value in the raster `afiel_path` is equal to value
        # calculated by hand
        result_set = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                afiel_path):
            result_set.update(numpy.unique(raster_block))
        assert len(result_set) == 1, (
            "One unique value expected in organic matter raster")
        test_result = list(result_set)[0]
        assert abs(test_result - 0.30895) < 0.0001, (
            "Test result does not match expected value")

    def test_calc_awilt(self):
        """Test the function _calc_awilt against value calculated by hand.

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

        def create_test_raster(target_path, fill_value):
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

        sand_path = os.path.join(self.workspace_dir, 'sand.tif')
        silt_path = os.path.join(self.workspace_dir, 'silt.tif')
        clay_path = os.path.join(self.workspace_dir, 'clay.tif')
        ompc_path = os.path.join(self.workspace_dir, 'ompc.tif')
        bulkd_path = os.path.join(self.workspace_dir, 'bulkd.tif')

        create_test_raster(sand_path, 0.39)
        create_test_raster(silt_path, 0.41)
        create_test_raster(clay_path, 0.2)
        create_test_raster(ompc_path, 0.913304)
        create_test_raster(bulkd_path, 1.5)

        awilt_path = os.path.join(self.workspace_dir, 'awilt.tif')

        forage._calc_awilt(
            sand_path, silt_path, clay_path, ompc_path, bulkd_path, awilt_path)

        # assert the value in the raster `awilt_path` is equal to value
        # calculated by hand
        result_set = set()
        for offset_map, raster_block in pygeoprocessing.iterblocks(
                awilt_path):
            result_set.update(numpy.unique(raster_block))
        assert len(result_set) == 1, (
            "One unique value expected in organic matter raster")
        test_result = list(result_set)[0]
        assert abs(test_result - 0.201988) < 0.0001, (
            "Test result does not match expected value")

    def test_afiel_awilt(self):
        """Test that values calculated by _afiel_awilt are reasonable.

        Use the function `_afiel_awilt` to calculate field capacity and wilting
        point from randomly generated inputs. Test that calculated field
        capacity and calculated wilting point are in the range [0.01, 0.9].
        Introduce nodata values into the input rasters and test that calculated
        field capacity and wilting point remain in the range [0.01, 0.9].

        Raises:
            AssertionError if calculated field capacity is outside the
                acceptable range
            AssertionError if calculated wilting point is outside the
                acceptable range

        Returns:
            None
        """
        from natcap.invest import forage

        def create_random_raster(target_path, lower_bound, upper_bound):
            """Create a small raster of random floats.

            The raster will have 10 rows and 10 columns and will be in the
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
            geotransform = [0, 1, 0, 44.5, 0, 1]
            n_cols = 10
            n_rows = 10
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

            random_array = numpy.random.uniform(
                lower_bound, upper_bound, (n_rows, n_cols))
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
                n_vals = numpy.random.randint(
                    0, (prior_copy.shape[0] * prior_copy.shape[1]))
                insertions = 0
                while insertions < n_vals:
                    row = numpy.random.randint(0, prior_copy.shape[0])
                    col = numpy.random.randint(0, prior_copy.shape[1])
                    modified_copy[row, col] = nodata_value
                    insertions += 1
                return modified_copy

            prior_copy = os.path.join(
                self.workspace_dir, 'prior_to_insert_nodata.tif')
            shutil.copyfile(target_raster, prior_copy)

            pygeoprocessing.raster_calculator(
                [(prior_copy, 1)],
                insert_op, target_raster,
                gdal.GDT_Float32, nodata_value)

            os.remove(prior_copy)

        def assert_all_values_in_raster_within_range(
                raster_to_test, minimum_acceptable_value,
                maximum_acceptable_value, nodata_value):
            """Test that raster contains values within acceptable range.

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
                assert min_val >= minimum_acceptable_value, (
                    "Raster contains values smaller than acceptable "
                    + "minimum: {}, {}".format(raster_to_test, min_val))
                max_val = numpy.amax(
                    raster_block[raster_block != nodata_value])
                assert max_val <= maximum_acceptable_value, (
                    "Raster contains values larger than acceptable "
                    + "maximum: {}, {}".format(raster_to_test, max_val))

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
        maximum_acceptable_value = 0.9
        nodata_value = _TARGET_NODATA

        forage._afiel_awilt(
            site_index_path, site_param_table, som1c_2_path,
            som2c_2_path, som3c_path, sand_path, silt_path, clay_path,
            bulk_d_path, pp_reg)

        for key, path in pp_reg.iteritems():
            assert_all_values_in_raster_within_range(
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
                assert_all_values_in_raster_within_range(
                    path, minimum_acceptable_value,
                    maximum_acceptable_value, nodata_value)

    def test_persistent_params(self):
        """Test that values calculated by persistent_params are reasonable.

        Use the function `persistent_params` to calculate wc, eftext, p1co2_2,
        fps1s3, and fps2s3 from randomly generated inputs. Test that each of
        the calculated quantities are within the range [0, 1].  Introduce
        nodata values into the inputs and test that calculated values
        remain inside the specified ranges.

        Raises:
            AssertionError if wc is outside the range [0, 1]
            AssertionError if eftext is outside the range [0, 1]
            AssertionError if p1co2_2 is outside the range [0, 1]
            AssertionError if fps1s3 is outside the range [0, 1]
            AssertionError if fps2s3 is outside the range [0, 1]

        Returns:
            None
        """
        from natcap.invest import forage

        def create_random_raster(target_path, lower_bound, upper_bound):
            """Create a small raster of random floats.

            The raster will have 10 rows and 10 columns and will be in the
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
            geotransform = [0, 1, 0, 44.5, 0, 1]
            n_cols = 10
            n_rows = 10
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

            random_array = numpy.random.uniform(
                lower_bound, upper_bound, (n_rows, n_cols))
            target_band.WriteArray(random_array)
            target_raster = None

        def insert_nodata_values_into_raster(target_raster, nodata_value):
            """Insert nodata at arbitrary locations in `target_raster`."""
            def insert_op(prior_copy):
                modified_copy = prior_copy
                n_vals = numpy.random.randint(
                    0, (prior_copy.shape[0] * prior_copy.shape[1]))
                insertions = 0
                while insertions < n_vals:
                    row = numpy.random.randint(0, prior_copy.shape[0])
                    col = numpy.random.randint(0, prior_copy.shape[1])
                    modified_copy[row, col] = nodata_value
                    insertions += 1
                return modified_copy

            prior_copy = os.path.join(
                self.workspace_dir, 'prior_to_insert_nodata.tif')
            shutil.copyfile(target_raster, prior_copy)

            pygeoprocessing.raster_calculator(
                [(prior_copy, 1)],
                insert_op, target_raster,
                gdal.GDT_Float32, nodata_value)

            os.remove(prior_copy)

        def assert_all_values_in_raster_within_range(
                raster_to_test, minimum_acceptable_value,
                maximum_acceptable_value, nodata_value):
            """Test that raster contains values within acceptable range.

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
                assert min_val >= minimum_acceptable_value, (
                    "Raster contains values smaller than acceptable "
                    + "minimum: {}, {}".format(raster_to_test, min_val))
                max_val = numpy.amax(
                    raster_block[raster_block != nodata_value])
                assert max_val <= maximum_acceptable_value, (
                    "Raster contains values larger than acceptable "
                    + "maximum: {}, {}".format(raster_to_test, max_val))

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
                'omlech_2': numpy.random.uniform(0.6, 0.18)},
                }

        pp_reg = {
            'afiel_1_path': os.path.join(self.workspace_dir, 'afiel_1.tif'),
            'awilt_1_path': os.path.join(self.workspace_dir, 'awilt.tif'),
            'wc_path': os.path.join(self.workspace_dir, 'wc.tif'),
            'eftext_path': os.path.join(self.workspace_dir, 'eftext.tif'),
            'p1co2_2_path': os.path.join(self.workspace_dir, 'p1co2_2.tif'),
            'fps1s3_path': os.path.join(self.workspace_dir, 'fps1s3.tif'),
            'fps2s3_path': os.path.join(self.workspace_dir, 'fps2s3.tif'),
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
                'minimum_acceptable_value': 0,
                'maximum_acceptable_value': 1,
                'nodata_value': _TARGET_NODATA,
                },
            'eftext_path': {
                'minimum_acceptable_value': 0,
                'maximum_acceptable_value': 1,
                'nodata_value': _IC_NODATA,
                },
            'p1co2_2_path': {
                'minimum_acceptable_value': 0,
                'maximum_acceptable_value': 1,
                'nodata_value': _IC_NODATA,
                },
            'fps1s3_path': {
                'minimum_acceptable_value': 0,
                'maximum_acceptable_value': 1,
                'nodata_value': _IC_NODATA,
                },
            'fps2s3_path': {
                'minimum_acceptable_value': 0,
                'maximum_acceptable_value': 1,
                'nodata_value': _IC_NODATA,
                },
        }

        forage._persistent_params(
            site_index_path, site_param_table, sand_path, clay_path, pp_reg)

        for path, ranges in acceptable_range_dict.iteritems():
            assert_all_values_in_raster_within_range(
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
                assert_all_values_in_raster_within_range(
                    pp_reg[path], ranges['minimum_acceptable_value'],
                    ranges['maximum_acceptable_value'],
                    ranges['nodata_value'])

    def test_aboveground_ratio(self):
        """Test that values calculated by `aboveground_ratio` are valid.

        Use the function `aboveground_ratio` to calculate the C/N or P
        ratio of decomposing aboveground material from random inputs. Test
        that the calculated ratio, agdrat, is within the range [1, 150].
        Introduce nodata values into the inputs and test that calculated
        agdrat remains inside the range [1, 150].

        Raises:
            AssertionError if agdrat is outside the range [1, 150]

        Returns:
            None
        """
        from natcap.invest import forage

        def assert_all_values_in_array_within_range(
                array_to_test, minimum_acceptable_value,
                maximum_acceptable_value, nodata_value):
            """Test that array contains values within acceptable range.

            The values within `array_to_test` that are not null must be
            greater than or equal to `minimum_acceptable_value` and
            less than or equal to `maximum_acceptable_value`.

            Raises:
                AssertionError if values are outside acceptable range

            Returns:
                None
            """
            if len(array_to_test[array_to_test != nodata_value]) == 0:
                pass
            min_val = numpy.amin(
                array_to_test[array_to_test != nodata_value])
            assert min_val >= minimum_acceptable_value, (
                "Array contains values smaller than acceptable minimum")
            max_val = numpy.amax(
                array_to_test[array_to_test != nodata_value])
            assert max_val <= maximum_acceptable_value, (
                "Array contains values larger than acceptable maximum")

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

        array_shape = (10, 10)

        tca = numpy.random.uniform(300, 700, array_shape)
        anps = numpy.random.uniform(1, numpy.amin(tca), array_shape)
        pcemic_1 = numpy.random.uniform(12, 20, array_shape)
        pcemic_2 = numpy.random.uniform(3, 11, array_shape)
        pcemic_3 = numpy.random.uniform(0.001, 0.1, array_shape)
        cemicb = (pcemic_2 - pcemic_1) / pcemic_3

        minimum_acceptable_agdrat = 0.01
        maximum_acceptable_agdrat = 1000
        agdrat_nodata = _TARGET_NODATA

        agdrat = forage._aboveground_ratio(
            anps, tca, pcemic_1, pcemic_2, pcemic_3, cemicb)

        assert_all_values_in_array_within_range(
            agdrat, minimum_acceptable_agdrat, maximum_acceptable_agdrat,
            agdrat_nodata)

        for input_array in [anps, tca]:
            insert_nodata_values_into_array(input_array, _TARGET_NODATA)
            agdrat = forage._aboveground_ratio(
                anps, tca, pcemic_1, pcemic_2, pcemic_3, cemicb)

            assert_all_values_in_array_within_range(
                agdrat, minimum_acceptable_agdrat, maximum_acceptable_agdrat,
                agdrat_nodata)
        for input_array in [pcemic_1, pcemic_2, pcemic_3, cemicb]:
            insert_nodata_values_into_array(input_array, _IC_NODATA)
            agdrat = forage._aboveground_ratio(
                anps, tca, pcemic_1, pcemic_2, pcemic_3, cemicb)

            assert_all_values_in_array_within_range(
                agdrat, minimum_acceptable_agdrat, maximum_acceptable_agdrat,
                agdrat_nodata)
