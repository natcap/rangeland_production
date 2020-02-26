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

A beta version of this model was described in Kowal et al. 2019,
available at https://www.nature.com/articles/s41598-019-56470-3.

For an interactive visualization of sample model outputs, see
http://viz.naturalcapitalproject.org/rangelands .

References
-------------------
Freer, M, A. D Moore, and J. R Donnelly. “The GRAZPLAN Animal Biology Model for Sheep and Cattle and the GrazFeed Decision Support Tool.” Canberra, ACT Australia: CSIRO Plant Industry, 2012.

Holland, E. A., Parton, W. J., Detling, J. K., and D.L. Coppock.  "Physiological Responses of Plant Populations to Herbivory and Their Consequences for Ecosystem Nutrient Flow." The American Naturalist 140, no. 4 (1992): 685-706. doi:10.1086/285435.

Kowal, V.A., Jones, S.M., Keesing, F., Allan, B., Schieltz, J., and R. Chaplin-Kramer. A coupled forage-grazer model predicts viability of livestock production and wildlife habitat at the regional scale. Scientific Repports 9, no. 19957 (2019). doi:10.1038/s41598-019-56470-3

Parton, W. J., J. M. O. Scurlock, D. S. Ojima, T. G. Gilmanov, R. J. Scholes, D. S. Schimel, T. Kirchner, et al. “Observations and Modeling of Biomass and Soil Organic Matter Dynamics for the Grassland Biome Worldwide.” Global Biogeochemical Cycles 7, no. 4 (1993): 785–809. doi:10.1029/93GB02042.

General Information
-------------------

* Website: https://naturalcapitalproject.org
* Source code: https://github.com/natcap/rangeland_production
* Issue tracker: https://github.com/natcap/rangeland_production/issues
* Users' guide: http://viz.naturalcapitalproject.org/rangelands/docs/Documentation_RPM_8.7.19.pdf

Dependencies
------------

Run ``make check`` to test if all required dependencies are installed on your system.


NSIS-specific requirements
++++++++++++++++++++++++++
The RPM NSIS installer requires the following:

* NSIS version < 3
* Installed Plugins:
    * Nsisunz: http://nsis.sourceforge.net/Nsisunz_plug-in
    * InetC: http://nsis.sourceforge.net/Inetc_plug-in
    * NsProcess: http://nsis.sourceforge.net/NsProcess_plugin

Managing python dependencies
++++++++++++++++++++++++++++
We recommend using a virtual environment to manage your python dependencies, and there is
a Makefile target to assist with this::

    $ make env
    $ source env/bin/activate

Or on Windows, use the following instead from a CMD prompt::

    > make env
    > .\env\bin\activate

This makefile target is included for convenience ... you may of course choose to
manage your own virtual environment.  ``requirements.txt``,
``requirements-dev.txt`` and ``requirements-gui.txt`` list the python
dependencies needed.

Using a different environment name
""""""""""""""""""""""""""""""""""
If you prefer a different name for your environment, you may pass the environment name as
a parameter to make::

    $ make ENV=myEnv env

You could then activate the environment created at ``myEnv``.


Using a different environment management tool
"""""""""""""""""""""""""""""""""""""""""""""
The RPM Makefile uses ``virtualenv`` to set up an environment, but this is
not the only `environment management tool out there
<https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments>`_.
You may elect to manage your virtual environment a different way, independent
of ``make env``.  The only requirement for the build process is that the required
tools are available on your PATH and the required python packages can be imported.


Building RPM Distributions
-----------------------------

Once the required tools and packages are available, we can build RPM.


Building ``rangeland_production`` python package
+++++++++++++++++++++++++++++++++++++++++

A Makefile target has been created for your convenience::

    $ make python_packages

This will create a wheel for your platform and a zip source archive in ``dist/``.
Both of these files (``dist/rangeland_production*.whl`` and ``dist/rangeland_production*.zip``)
can be installed by pip.

Building python packages without GNU make
"""""""""""""""""""""""""""""""""""""""""
Python distributions may be built with the standard distutils/setuptools commands::

    $ python setup.py bdist_wheel
    $ python setup.py sdist

RPM Standalone Binaries
++++++++++++++++++++++++++

Once the appropriate dependencies are available, RPM can also be built as a
standalone application::

    $ make binaries

An important detail about building binaries is that ``rangeland_production`` must be
installed as a wheel to ensure that the distribution information is in the
correct location.

This will create a directory at ``dist/rangeland_production`` holding the application
binaries and relevant shared libraries.

Binaries cannot be cross-compiled for other operating systems.


RPM Windows Installer
++++++++++++++++++++++++

The RPM installer for Windows can be built with::

    > make windows_installer

This will create the installer at ``dist/rangeland_production*_Setup.exe``.


Tests
-----

RPM includes a suite of unit tests to ensure software quality.

Model tests
+++++++++++

To run tests on RPM::

    $ make test


Changing how GNU make runs tests
++++++++++++++++++++++++++++++++

The InVEST Makefile setup depends on ``nosetests`` and takes advantage of its
plugins for line coverage and xunit reports.  You can force ``make`` to use a
different test runner by setting a parameter at the command line.  For example,
to run the tests with ``pytest``::

    $ make TESTRUNNER=pytest test


Copyright and license information
---------------------------------

A file called ``LICENSE.txt`` should have accompanied this distribution.  If it
is missing, the license may be found on our project page,
https://github.com/natcap/rangeland_production
