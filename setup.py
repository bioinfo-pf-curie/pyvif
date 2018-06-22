#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

_MAJOR               = 0
_MINOR               = 1
_MICRO               = 0
version              = '{}.{}.{}'.format(_MAJOR, _MINOR, _MICRO)
release              = '{}.{}'.format(_MAJOR, _MINOR)

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=6.0', 'jinja2', 'matplotlib', 'numpy', 'pandas',
                'pysam', 'pytest', ]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    name='pyvif',
    version=version,
    author="Dimitri Desvillechabrol",
    author_email='dimitri.desvillechabrol@curie.fr',
    maintainer="Dimitri Desvillechabrol",
    maintainer_email='dimitri.desvillechabrol@curie.fr',
    keywords=['pyvif', 'pacbio', 'virus', 'breakpoints'],
    description="Python Virus Integration Finder detects integration site of"
                " virus in human genome using capture pacbio.",
    long_description=readme + '\n\n' + history,
    license="BSD license",
    plateforms=['Linux', 'Unix', 'MacOsX'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=requirements,
    include_package_data=True,
    packages=find_packages(include=['pyvif']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://bitbucket.org/ddesvillechabrol/pyvif',
    zip_safe=False,
    entry_points={
        'console_scripts': [
        ],
    },
)
