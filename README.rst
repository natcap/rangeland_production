RPM: Rangeland Production Model
================================================================

The Rangeland Production Model is a simulation that represents dynamic growth
of grass and consumption of grass by herbivores to predict productivity in
terms of herbivore diet sufficiency, or the predicted intake of protein and
energy relative to maintenance requirements. The model incorporates a gridded
(i.e., pixel- or raster-based) implementation of the Century ecosystem model
(Parton et al. 1993) coupled with a basic physiology submodel adapted from
GRAZPLAN (Freer et al. 2012). At each monthly timestep for which it is run,
the RPM uses the gridded implementation of Century to predict forage biomass
and nutrient content on each pixel covering the study area. The model
estimates the pixel-level densities of herbivores, estimates diet selection
of those herbivores, and integrates the removal of forage by herbivores into
the production of forage biomass in the following timestep.

For full model documentation and a set of sample outputs, see 
http://viz.naturalcapitalproject.org/rangelands .

General Information
-------------------

* Website: https://naturalcapitalproject.org
* Source code: https://github.com/natcap/rangeland_production
* Issue tracker: https://github.com/natcap/rangeland_production/issues
* Users' guide: http://viz.naturalcapitalproject.org/rangelands/docs/Documentation_RPM_8.7.19.pdf

Dependencies
------------

Run ``make check`` to test if all required dependencies are installed on your system.
OS-specific installation instructions are found either online at
http://invest.readthedocs.io/en/latest/installing.html or locally at ``doc/api-docs/installing.rst``.


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
