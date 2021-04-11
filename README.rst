RPM: Rangeland Production Model
================================================================

The purpose of the Rangeland Production model is to project forage
production and diet sufficiency for grazing animals under different
conditions of climate and management. The model incorporates a gridded
(i.e., pixel- or raster-based) implementation of the Century ecosystem
model (Parton et al. 1993) coupled with a basic physiology submodel
adapted from GRAZPLAN (Freer et al. 2012). The raster-based Century
implementation submodel simulates the growth of herbaceous forage on
each pixel of the simulated area, according to climate and soil
conditions, at a monthly timestep. The ruminant physiology submodel
adapted from GRAZPLAN calculates offtake of forage by grazing animals
according to the biomass and protein content of the simulated forage,
and estimates the adequacy of the diet to meet the animals’ energy
requirements.  Then, the estimated offtake by animals is integrated
into the regrowth of forage in the following timestep through impacts
on simulated potential production, root:shoot ratio, and plant nitrogen
content according to Century’s existing grazing routine (Holland et al.
1992).

Inputs to the Rangeland model therefore include gridded inputs to
Century and descriptive parameters of the modeled livestock herd.
Outputs of the model consist of monthly rasters giving forage biomass and
protein content, forage intake by grazing animals, and an estimate of the
adequacy of the forage consumed by animals to meet their maintenance energy
requirements, over the same time period as the inputs.

The model was described by Kowal et al. 2021, available at 
https://doi.org/10.3390/land10040397. An earlier, point-based version of
the model was described in Kowal et al. 2019, available at
https://www.nature.com/articles/s41598-019-56470-3.

For an interactive visualization of sample model outputs, see
http://viz.naturalcapitalproject.org/rangelands .

References
-------------------
Freer, M, A. D Moore, and J. R Donnelly. “The GRAZPLAN Animal Biology Model for Sheep and Cattle and the GrazFeed Decision Support Tool.” Canberra, ACT Australia: CSIRO Plant Industry, 2012.

Holland, E. A., Parton, W. J., Detling, J. K., and D.L. Coppock.  "Physiological Responses of Plant Populations to Herbivory and Their Consequences for Ecosystem Nutrient Flow." The American Naturalist 140, no. 4 (1992): 685-706. doi:10.1086/285435.

Kowal, V.A., Jones, S.M., Keesing, F., Allan, B., Schieltz, J., and R. Chaplin-Kramer. "A coupled forage-grazer model predicts viability of livestock production and wildlife habitat at the regional scale." Scientific Repports 9, no. 19957 (2019). doi:10.1038/s41598-019-56470-3

Kowal, V.A., Ahlborn, J., Jamsranjav, C., Avirmed, O., and R. Chaplin-Kramer. "Modeling integrated impacts of climate change and grazing on Mongolia’s rangelands." Land 10, no. 397 (2021). doi:10.3390/land10040397

Parton, W. J., J. M. O. Scurlock, D. S. Ojima, T. G. Gilmanov, R. J. Scholes, D. S. Schimel, T. Kirchner, et al. “Observations and Modeling of Biomass and Soil Organic Matter Dynamics for the Grassland Biome Worldwide.” Global Biogeochemical Cycles 7, no. 4 (1993): 785–809. doi:10.1029/93GB02042.

General Information
-------------------

* Installer, sample data, and complete users' guide: https://github.com/natcap/rangeland_production/releases/latest
* Website: https://naturalcapitalproject.org
* Source code: https://github.com/natcap/rangeland_production
* Issue tracker: https://github.com/natcap/rangeland_production/issues

Copyright and license information
---------------------------------

A file called ``LICENSE.txt`` should have accompanied this distribution.  If it
is missing, the license may be found on our project page,
https://github.com/natcap/rangeland_production
