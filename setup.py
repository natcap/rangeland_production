"""setup.py module for rangeland_production

Common functionality provided by setup.py:
    build_sphinx

For other commands, try `python setup.py --help-commands`
"""
from setuptools.extension import Extension
from setuptools import setup


# Read in requirements.txt and populate the python readme with the
# non-comment, non-environment-specifier contents.
_REQUIREMENTS = [req.split(';')[0].split('#')[0].strip() for req in
                 open('requirements.txt').readlines()
                 if not req.startswith(('#', 'hg+')) and len(req.strip()) > 0]
_GUI_REQUIREMENTS = [req.split(';')[0].split('#')[0].strip() for req in
                     open('requirements-gui.txt').readlines()
                     if not req.startswith(('#', 'hg+')) and len(req.strip()) > 0]
README = open('README_PYTHON.rst').read().format(
    requirements='\n'.join(['    ' + r for r in _REQUIREMENTS]))


setup(
    name='rangeland_production',
    description="Rangeland Production Model",
    long_description=README,
    maintainer='Ginger Kowal',
    maintainer_email='gkowal@stanford.edu',
    url='http://bitbucket.org/natcap/rangeland_production',
    packages=[
        'rangeland_production',
        'rangeland_production.ui',
    ],
    package_dir={
        'rangeland_production': 'src/rangeland_production'
    },
    use_scm_version={'version_scheme': 'post-release',
                     'local_scheme': 'node-and-date'},
    include_package_data=True,
    install_requires=_REQUIREMENTS,
    setup_requires=['setuptools_scm'],
    license='BSD',
    zip_safe=False,
    keywords='gis rangeland',
    classifiers=[
        'Intended Audience :: Developers',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: Microsoft',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: GIS'
    ],
    entry_points={
        'console_scripts': [
            'rangeland_production = rangeland_production.cli:main'
        ],
    },
    extras_require={
        'ui': _GUI_REQUIREMENTS,
    },
)
