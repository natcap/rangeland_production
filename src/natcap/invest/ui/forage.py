# coding=UTF-8

from natcap.invest.ui import model, inputs
import natcap.invest.forage


class forage(model.InVESTModel):
    def __init__(self):
        model.InVESTModel.__init__(
            self,
            label=u'Rangeland Production',
            target=forage.execute,
            validator=forage.validate,
            localdoc=u'forage.html'),

        self.n_months = inputs.Text(
            args_key=u'n_months',
            helptext=(
                u"Number of months for which the model should be run "
                "(1 or greater)."),
            label=u'Number of Months',
            validator=self.validator)
        self.add_input(self.n_months)
        self.starting_year = inputs.Text(
            args_key=u'starting_year',
            helptext=(
                u"The first year of the model run.  For example, if the first "
                "month of the simulation is January of 2010, this value "
                "should be 2010."),
            label=u'Starting Year',
            validator=self.validator)
        self.add_input(self.starting_year)
        self.starting_month = inputs.Text(
            args_key=u'starting_month',
            helptext=(
                u"The first month of the model run (1-12). For example, if "
                "the first month of the model run is January, this value "
                "should be 1."),
            label=u'Starting month',
            validator=self.validator)
        self.add_input(self.starting_month)
        self.aoi_path = inputs.File(
            args_key=u'aoi_path',
            helptext=(
                u"An OGR-supported, single-feature polygon vector in "
                "geographic coordinates (i.e., latitude and longitude). All "
                "model outputs will be clipped to the extent of the feature "
                "in this vector."),
            label=u'Area of Interest (Vector)',
            validator=self.validator)
        self.add_input(self.aoi_path)
        self.management_threshold = inputs.Text(
            args_key=u'management_threshold',
            helptext=(
                u"Total biomass (kg/ha), including live and standing dead "
                "biomass across all plant functional types, that must be left "
                "standing after offtake by grazing animals. Total offtake by "
                "animals is limited at each model step such that residual "
                "biomass is equal to or greater than this value."),
            label=u'Management Threshold',
            validator=self.validator)
        self.add_input(self.management_threshold)
        self.proportion_legume_path = inputs.File(
            args_key=u'proportion_legume_path',
            helptext=(
                u"A GDAL-supported raster file containing the proportion of "
                "legume of dry aboveground biomass, by weight (0-1). For "
                "example, if 100% of pasture biomass is constituted by "
                "legumes, this value should be 1.  Must be in geographic "
                "coordinates (latitude and longitude)."),
            label=u'Proportion Legume (Raster)',
            validator=self.validator)
        self.add_input(self.proportion_legume_path)
        self.clay_proportion_path = inputs.File(
            args_key=u'clay_proportion_path',
            helptext=(
                u"A GDAL-supported raster file giving the fraction of soil "
                "that is clay (0-1). Must be in geographic coordinates "
                "(latitude and longitude)."),
            label=u'Proportion Clay (Raster)',
            validator=self.validator)
        self.add_input(self.clay_proportion_path)
        self.silt_proportion_path = inputs.File(
            args_key=u'silt_proportion_path',
            helptext=(
                u"A GDAL-supported raster file giving the fraction of soil "
                "that is silt (0-1). Must be in geographic coordinates "
                "(latitude and longitude)."),
            label=u'Proportion Silt (Raster)',
            validator=self.validator)
        self.add_input(self.silt_proportion_path)
        self.sand_proportion_path = inputs.File(
            args_key=u'sand_proportion_path',
            helptext=(
                u"A GDAL-supported raster file giving the fraction of soil "
                "that is sand (0-1). Must be in geographic coordinates "
                "(latitude and longitude)."),
            label=u'Proportion Sand (Raster)',
            validator=self.validator)
        self.add_input(self.sand_proportion_path)
        self.bulk_density_path = inputs.File(
            args_key=u'bulk_density_path',
            helptext=(
                u"A GDAL-supported raster file containing the bulk density of "
                "soil, in grams per cubic centimeter. Must be in geographic "
                "coordinates (latitude and longitude)."),
            label=u'Bulk Density (Raster)',
            validator=self.validator)
        self.add_input(self.bulk_density_path)
        self.ph_path = inputs.File(
            args_key=u'ph_path',
            helptext=(
                u"A GDAL-supported raster file containing soil pH (0-14). "
                "Must be in geographic coordinates (latitude and longitude)."),
            label=u'Soil pH (Raster)',
            validator=self.validator)
        self.add_input(self.ph_path)
        self.min_temp_path_pattern = inputs.Text(
            args_key=u'min_temp_path_pattern',
            helptext=(
                u"An example file path to monthly average daily minimum "
                "temperature (degrees C), which must contain the substring "
                "'<month>'. The model only requires one raster of minimum "
                "temperature per month, even if the model is run for multiple "
                "years.  For example, if this value is given as "
                "'./workspace/min_temp/<month>.tif', Starting Year is 2010, "
                "Starting Month is 1, and Number of Months is 13, the model "
                "will expect to find raster files named "
                "'./workspace/min_temp/01.tif' to "
                "'./workspace/min_temp/12.tif'.  The '<month>' value in "
                "input files must be two digits.  The rasters must be in "
                "geographic coordinates."),
            label=u'Monthly Minimum Temperature Path Pattern',
            validator=self.validator)
        self.add_input(self.min_temp_path_pattern)
        self.max_temp_path_pattern = inputs.Text(
            args_key=u'max_temp_path_pattern',
            helptext=(
                u"An example file path to monthly average daily maximum "
                "temperature (degrees C), which must contain the substring "
                "'<month>'. The model only requires one raster of maximum "
                "temperature per month, even if the model is run for multiple "
                "years.  For example, if this value is given as "
                "'./workspace/max_temp/<month>.tif', Starting Year is 2010, "
                "Starting Month is 1, and Number of Months is 13, the model "
                "will expect to find raster files named "
                "'./workspace/max_temp/01.tif' to "
                "'./workspace/max_temp/12.tif'.  The '<month>' value in "
                "input files must be two digits.  The rasters must be in "
                "geographic coordinates."),
            label=u'Monthly Maximum Temperature Path Pattern',
            validator=self.validator)
        self.add_input(self.max_temp_path_pattern)
        self.monthly_precip_path_pattern = inputs.Text(
            args_key=u'monthly_precip_path_pattern',
            helptext=(
                u"An example file path to monthly total precipitation rasters "
                "(cm), which must contain the substrings '<month>' and "
                "'<year>'. At least 12 months of precipitation must be "
                "supplied, and precipitation must be supplied for each month "
                "that the model is run. For example, if this value is given "
                "as './workspace/precip/<month>_<year>.tif', Starting Year "
                "is 2010, Starting Month is 1, and Number of Months is 13, "
                "the model will expect to find raster files named "
                "'./workspace/precip/01_2010.tif' to "
                "'./workspace/precip/01_2011.tif'.  The '<month>' value "
                "in input files must be two digits.  The rasters must be in "
                "geographic coordinates."),
            label=u'Monthly Precipitation Path Pattern',
            validator=self.validator)
        self.add_input(self.monthly_precip_path_pattern)
        self.monthly_vi_path_pattern = inputs.Text(
            args_key=u'monthly_vi_path_pattern',
            helptext=(
                u"An example file path to monthly vegetation index rasters "
                "such as NDVI, which must contain the substrings '<month>' "
                "and '<year>'.  One vegetation index raster must be "
                "supplied for each month that the model is run. For example, "
                "if this value is given as "
                "'./workspace/NDVI/<month>_<year>.tif', Starting Year is "
                "2010, Starting Month is 1, and Number of Months is 13, the "
                "model will expect to find raster files named "
                "'./workspace/NDVI/01_2010.tif' to "
                "'./workspace/NDVI/01_2011.tif'.  The '<month>' value in "
                "input files must be two digits.  The rasters must be in "
                "geographic coordinates."),
            label=u'Monthly Vegetation Index Path Pattern',
            validator=self.validator)
        self.add_input(self.monthly_vi_path_pattern)
        self.site_param_spatial_index_path = inputs.File(
            args_key=u'site_param_spatial_index_path',
            helptext=(
                u"A GDAL-supported raster of integers giving site codes for "
                "the area of interest. These site codes must match values in "
                "the "site" column of the Site Parameter Table. The raster "
                "must contain integer values.  The raster must be in "
                "geographic coordinates."),
            label=u'Site Spatial Index (Raster)',
            validator=self.validator)
        self.add_input(self.site_param_spatial_index_path)
        self.veg_spatial_composition_path_pattern = inputs.Text(
            args_key=u'veg_spatial_composition_path_pattern',
            helptext=(
                u"An example file path to vegetation fractional cover "
                "rasters, one per plant functional type.  The path must "
                "contain the substring '<PFT>', where '<PFT>' can be "
                "replaced with an integer value for each simulated plant "
                "functional type. Integer values for each plant functional "
                "type must match values in the column named 'PFT' in the "
                "plant functional type parameter table.  For example, if this "
                "value is given as './workspace/pft_<PFT>.tif' and the "
                "directory 'workspace' contains these files: "
                "'pft_1.tif', 'pft_10.tif', then the 'PFT' column of "
                "the plant functional type parameter table must contain the "
                "values 1 and 10.  Fractional cover rasters must be in "
                "geographic coordinates."),
            label=u'Plant Functional Type Fractional Cover Pattern',
            validator=self.validator)
        self.add_input(self.veg_spatial_composition_path_pattern)
        self.animal_grazing_areas_path = inputs.File(
            args_key=u'animal_grazing_areas_path',
            helptext=(
                u"An OGR-supported polygon vector in geographic coordinates "
                "(i.e., latitude and longitude), giving the areas over which "
                "animals graze.  Each polygon feature may contain one animal "
                "type, and features may not overlap. The vector must contain "
                "a column named 'animal_id' with integer values that match "
                "values given in the 'animal_id' column of the animal "
                "parameter table, and a column named 'num_animal' that "
                "gives the number of animals grazing inside each polygon "
                "feature."),
            label=u'Animal Grazing Areas (Vector)',
            validator=self.validator)
        self.add_input(self.animal_grazing_areas_path)
        self.site_param_table = inputs.File(
            args_key=u'site_param_table',
            helptext=(
                u"A CSV table giving parameters for each site code indexed by "
                "the site spatial index raster. The table must contain a "
                "column named 'site' with integer values matching values in "
                "the site spatial index raster, thereby linking parameters of "
                "each modeled site type to their location within the study "
                "area. See model documentation for other required parameters "
                "that must be included in this table."),
            label=u'Site Parameter Table (CSV)',
            validator=self.validator)
        self.add_input(self.site_param_table)
        self.veg_trait_path = inputs.File(
            args_key=u'veg_trait_path',
            helptext=(
                u"A CSV table giving parameters for each plant functional t"
                "ype. The table must contain a column named 'PFT' with "
                "integer values that match integer values in corresponding "
                "plant functional type fractional cover raster file names, "
                "thereby linking plant functional type parameters to their "
                "fractional cover in each pixel. See model documentation for "
                "other required parameters that must be included in this "
                "table."),
            label=u'Plant Functional Type Parameter Table (CSV)',
            validator=self.validator)
        self.add_input(self.veg_trait_path)
        self.animal_trait_path = inputs.File(
            args_key=u'animal_trait_path',
            helptext=(
                u"A CSV table giving parameters for each grazing animal type. "
                "The table must contain a column named 'animal_id' "
                "containing integer values that match integer values in the "
                "'animal_id' field of the animal grazing areas vector "
                "layer, thereby linking parameters of each animal type to the "
                "area over which they graze. See model documentation for "
                "other required parameters that must be included in this "
                "table."),
            label=u'Animal Parameter Table (CSV)',
            validator=self.validator)
        self.add_input(self.animal_trait_path)
        self.site_initial_table = inputs.File(
            args_key=u'site_initial_table',
            helptext=(
                u"A CSV table giving initial values for each site-level state "
                "variable. The table must contain a column named 'site' "
                "with integer values matching values in the site spatial "
                "index raster, thereby linking initial values for each "
                "modeled site type to their location within the study area. "
                "See model documentation for state variables that are "
                "required to be included in this table."),
            label=u'Initial Conditions Site State Variables',
            validator=self.validator)
        self.add_input(self.site_initial_table)
        self.pft_initial_table = inputs.File(
            args_key=u'pft_initial_table',
            helptext=(
                u"A CSV table giving initial values for each PFT-level state "
                "variable. The table must contain a column named 'PFT' with "
                "integer values that match integer values in corresponding "
                "plant functional type fractional cover raster file names, "
                "thereby linking plant functional type initial conditions to "
                "their fractional cover in each pixel. See model "
                "documentation for state variables that are required to be "
                "included in this table."),
            label=u'Initial Conditions: Plant Functional Type State Variables',
            validator=self.validator)
        self.add_input(self.pft_initial_table)
        self.results_suffix = inputs.Text(
            args_key=u'results_suffix',
            helptext=(
                u"A string that will be added to the end of the model output "
                "files."),
            label=u'Results Suffix',
            validator=self.validator)
        self.add_input(self.results_suffix)

    def assemble_args(self):
        args = {
            self.workspace_dir.args_key: self.workspace_dir.value(),
            self.n_months.args_key: self.n_months.value(),
            self.starting_year.args_key: self.starting_year.value(),
            self.starting_month.args_key: self.starting_month.value(),
            self.aoi_path.args_key: self.aoi_path.value(),
            self.management_threshold.args_key:
                self.management_threshold.value(),
            self.proportion_legume_path.args_key:
                self.proportion_legume_path.value(),
            self.clay_proportion_path.args_key:
                self.clay_proportion_path.value(),
            self.silt_proportion_path.args_key:
                self.silt_proportion_path.value(),
            self.sand_proportion_path.args_key:
                self.sand_proportion_path.value(),
            self.bulk_density_path.args_key: self.bulk_density_path.value(),
            self.ph_path.args_key: self.ph_path.value(),
            self.min_temp_path_pattern.args_key:
                self.min_temp_path_pattern.value(),
            self.max_temp_path_pattern.args_key:
                self.max_temp_path_pattern.value(),
            self.monthly_precip_path_pattern.args_key:
                self.monthly_precip_path_pattern.value(),
            self..args_key: self..value(),
            self.site_param_spatial_index_path.args_key:
                self.site_param_spatial_index_path.value(),
            self.veg_spatial_composition_path_pattern.args_key:
                self.veg_spatial_composition_path_pattern.value(),
            self.animal_grazing_areas_path.args_key:
                self.animal_grazing_areas_path.value(),
            self.site_param_table.args_key: self.site_param_table.value(),
            self.veg_trait_path.args_key: self.veg_trait_path.value(),
            self.animal_trait_path.args_key: self.animal_trait_path.value(),
            self.site_initial_table.args_key: self.site_initial_table.value(),
            self.pft_initial_table.args_key: self.pft_initial_table.value(),
            self.results_suffix.args_key: self.results_suffix.value(),
        }

        return args
