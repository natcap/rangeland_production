# coding=UTF-8
import functools

from natcap.invest.ui import model, inputs
from natcap.invest.coastal_blue_carbon import coastal_blue_carbon
from natcap.invest.coastal_blue_carbon import preprocessor


class CoastalBlueCarbonPreprocessor(model.Model):
    def __init__(self):
        model.Model.__init__(
            self,
            label=u'Coastal Blue Carbon Preprocessor',
            target=preprocessor.execute,
            validator=preprocessor.validate,
            localdoc=u'../documentation/coastal_blue_carbon.html')

        self.lulc_lookup_uri = inputs.File(
            args_key=u'lulc_lookup_uri',
            helptext=(
                u"A CSV table used to map lulc classes to their values "
                u"in a raster, as well as to indicate whether or not "
                u"the lulc class is a coastal blue carbon habitat."),
            label=u'LULC Lookup Table (CSV)',
            validator=self.validator)
        self.add_input(self.lulc_lookup_uri)
        self.lulc_snapshot_list = inputs.Multi(
            args_key=u'lulc_snapshot_list',
            callable_=functools.partial(inputs.File, label="Input"),
            label=u'Land Use/Land Cover Rasters (GDAL-supported)',
            link_text=u'Add Another')
        self.add_input(self.lulc_snapshot_list)

    def assemble_args(self):
        args = {
            self.workspace.args_key: self.workspace.value(),
            self.suffix.args_key: self.suffix.value(),
            self.lulc_lookup_uri.args_key: self.lulc_lookup_uri.value(),
            self.lulc_snapshot_list.args_key: self.lulc_snapshot_list.value(),
        }
        return args


class CoastalBlueCarbon(model.Model):
    def __init__(self):
        model.Model.__init__(
            self,
            label=u'Coastal Blue Carbon',
            target=coastal_blue_carbon.execute,
            validator=coastal_blue_carbon.validate,
            localdoc=u'../documentation/coastal_blue_carbon.html')

        self.lulc_lookup_uri = inputs.File(
            args_key=u'lulc_lookup_uri',
            helptext=(
                u"A CSV table used to map lulc classes to their values "
                u"in a raster and to indicate whether or not the lulc "
                u"class is a coastal blue carbon habitat."),
            label=u'LULC Lookup Table (CSV)',
            validator=self.validator)
        self.add_input(self.lulc_lookup_uri)
        self.lulc_transition_matrix_uri = inputs.File(
            args_key=u'lulc_transition_matrix_uri',
            helptext=(
                u"Generated by the preprocessor.  This file must be "
                u"edited before it can be used by the main model.  The "
                u"left-most column represents the source lulc class, "
                u"and the top row represents the destination lulc "
                u"class."),
            label=u'LULC Transition Effect of Carbon Table (CSV)',
            validator=self.validator)
        self.add_input(self.lulc_transition_matrix_uri)
        self.carbon_pool_initial_uri = inputs.File(
            args_key=u'carbon_pool_initial_uri',
            helptext=(
                u"The provided CSV table contains information related "
                u"to the initial conditions of the carbon stock within "
                u"each of the three pools of a habitat.  Biomass "
                u"includes carbon stored above and below ground.  All "
                u"non-coastal blue carbon habitat lulc classes are "
                u"assumed to contain no carbon.  The values for "
                u"‘biomass’, ‘soil’, and ‘litter’ should be given in "
                u"terms of Megatonnes CO2e/ha."),
            label=u'Carbon Pool Initial Variables Table (CSV)',
            validator=self.validator)
        self.add_input(self.carbon_pool_initial_uri)
        self.carbon_pool_transient_uri = inputs.File(
            args_key=u'carbon_pool_transient_uri',
            helptext=(
                u"The provided CSV table contains information related "
                u"to the transition of carbon into and out of coastal "
                u"blue carbon pools.  All non-coastal blue carbon "
                u"habitat lulc classes are assumed to neither sequester "
                u"nor emit carbon as a result of change.  The "
                u"‘yearly_accumulation’ values should be given in terms "
                u"of Megatonnes of CO2e/ha-yr.  The ‘half-life’ values "
                u"must be given in terms of years.  The ‘disturbance’ "
                u"values must be given as a decimal percentage of stock "
                u"distrubed given a transition occurs away from a lulc- "
                u"class."),
            label=u'Carbon Pool Transient Variables Table (CSV)',
            validator=self.validator)
        self.add_input(self.carbon_pool_transient_uri)
        self.lulc_baseline_map_uri = inputs.File(
            args_key=u'lulc_baseline_map_uri',
            helptext=(
                u"A GDAL-supported raster representing the baseline "
                u"landscape/seascape."),
            label=u'Baseline LULC Raster (GDAL-supported)',
            validator=self.validator)
        self.add_input(self.lulc_baseline_map_uri)
        self.lulc_baseline_year = inputs.Text(
            args_key=u'lulc_baseline_year',
            label=u'Year of baseline LULC raster',
            validator=self.validator)
        self.add_input(self.lulc_baseline_year)
        self.lulc_transition_maps_list = inputs.Multi(
            args_key=u'lulc_transition_maps_list',
            callable_=functools.partial(inputs.File, label="Input"),
            label=u'Transition LULC Rasters (GDAL-supported)',
            link_text=u'Add Another')
        self.add_input(self.lulc_transition_maps_list)
        self.lulc_transition_years_list = inputs.Multi(
            args_key=u'lulc_transition_years_list',
            callable_=functools.partial(inputs.Text, label="Input"),
            label=u'Transition Years',
            link_text=u'Add Another')
        self.add_input(self.lulc_transition_years_list)
        self.analysis_year = inputs.Text(
            args_key=u'analysis_year',
            helptext=(
                u"An analysis year extends the transient analysis "
                u"beyond the transition years."),
            label=u'Analysis Year (Optional)',
            validator=self.validator)
        self.add_input(self.analysis_year)
        self.do_economic_analysis = inputs.Container(
            args_key=u'do_economic_analysis',
            expandable=True,
            expanded=True,
            label=u'Calculate Net Present Value of Sequestered Carbon')
        self.add_input(self.do_economic_analysis)
        self.do_price_table = inputs.Checkbox(
            args_key=u'do_price_table',
            helptext=u'',
            label=u'Use Price Table')
        self.do_economic_analysis.add_input(self.do_price_table)
        self.price = inputs.Text(
            args_key=u'price',
            helptext=u'The price per Megatonne CO2e at the base year.',
            label=u'Price',
            validator=self.validator)
        self.do_economic_analysis.add_input(self.price)
        self.interest_rate = inputs.Text(
            args_key=u'interest_rate',
            helptext=(
                u"The interest rate on the price per Megatonne CO2e, "
                u"compounded yearly."),
            label=u'Interest Rate (%)',
            validator=self.validator)
        self.do_economic_analysis.add_input(self.interest_rate)
        self.price_table_uri = inputs.File(
            args_key=u'price_table_uri',
            helptext=(
                u"Can be used in place of price and interest rate "
                u"inputs.  The provided CSV table contains the price "
                u"per Megatonne CO2e sequestered for a given year, for "
                u"all years from the original snapshot to the analysis "
                u"year, if provided."),
            interactive=False,
            label=u'Price Table (CSV)',
            validator=self.validator)
        self.do_economic_analysis.add_input(self.price_table_uri)
        self.discount_rate = inputs.Text(
            args_key=u'discount_rate',
            helptext=(
                u"The discount rate on future valuations of "
                u"sequestered carbon, compounded yearly."),
            label=u'Discount Rate (%)',
            validator=self.validator)
        self.do_economic_analysis.add_input(self.discount_rate)

        # Set interactivity, requirement as input sufficiency changes
        self.do_price_table.sufficiency_changed.connect(
            self.price.set_noninteractive)
        self.do_price_table.sufficiency_changed.connect(
            self.interest_rate.set_noninteractive)
        self.do_price_table.sufficiency_changed.connect(
            self.price_table_uri.set_interactive)

    def assemble_args(self):
        args = {
            self.workspace.args_key: self.workspace.value(),
            self.suffix.args_key: self.suffix.value(),
            self.lulc_lookup_uri.args_key: self.lulc_lookup_uri.value(),
            self.lulc_transition_matrix_uri.args_key:
                self.lulc_transition_matrix_uri.value(),
            self.carbon_pool_initial_uri.args_key:
                self.carbon_pool_initial_uri.value(),
            self.carbon_pool_transient_uri.args_key:
                self.carbon_pool_transient_uri.value(),
            self.lulc_baseline_map_uri.args_key:
                self.lulc_baseline_map_uri.value(),
            self.lulc_baseline_year.args_key:
                self.lulc_baseline_year.value(),
            self.lulc_transition_maps_list.args_key:
                self.lulc_transition_maps_list.value(),
            self.lulc_transition_years_list.args_key:
                self.lulc_transition_years_list.value(),
            self.analysis_year.args_key: self.analysis_year.value(),
            self.do_economic_analysis.args_key:
                self.do_economic_analysis.value(),
        }

        if self.do_economic_analysis.value():
            args[self.do_price_table.args_key] = self.do_price_table.value()
            args[self.price.args_key] = self.price.value()
            args[self.interest_rate.args_key] = self.interest_rate.value()
            args[self.price_table_uri.args_key] = self.price_table_uri.value()
            args[self.discount_rate.args_key] = self.discount_rate.value()

        return args
