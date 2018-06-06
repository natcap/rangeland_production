"""InVEST forage model tests."""

import unittest
import tempfile
import shutil
import os
import sys


SAMPLE_DATA = "C:/Users/ginge/Documents/NatCap/sample_inputs"
REGRESSION_DATA = "C:/Users/ginge/Documents/NatCap/regression_test_data"

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
            'bulk_density_path': os.path.join(SAMPLE_DATA,
                'soil', 'bulkd.tif'),
            'ph_path': os.path.join(SAMPLE_DATA,
                'soil', 'phihox_sl3.tif'),
            'clay_proportion_path': os.path.join(SAMPLE_DATA,
                'soil', 'clay.tif'),
            'silt_proportion_path': os.path.join(SAMPLE_DATA,
                'soil', 'silt.tif'),
            'sand_proportion_path': os.path.join(SAMPLE_DATA,
                'soil', 'sand.tif'),
            'monthly_precip_path_pattern': os.path.join(
                SAMPLE_DATA, 'CHIRPS_div_by_10',
                'chirps-v2.0.<year>.<month>.tif'),
            'min_temp_path_pattern': os.path.join(
                SAMPLE_DATA, 'temp', 'wc2.0_30s_tmin_<month>.tif'),
            'max_temp_path_pattern': os.path.join(
                SAMPLE_DATA, 'temp', 'wc2.0_30s_tmax_<month>.tif'),
            'site_param_table': os.path.join(SAMPLE_DATA,
                'site_parameters.csv'),
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
        """Launch forage model, ensure it runs!"""
        from natcap.invest import forage
        
        if not os.path.exists(SAMPLE_DATA):
            self.fail(
                "Sample input directory not found at %s" % SAMPLE_DATA)
                
        args = foragetests.generate_base_args(self.workspace_dir)
        forage.execute(args)
