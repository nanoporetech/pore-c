#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import versioneer

from setuptools import find_packages
from setuptools import setup

version = versioneer.get_version()

python_requires = ">=3.6"

install_requires = [
	"biopython ",
	"click >=7.0",
	"cooler >=0.8.5",
	"dask >=2.0.0",
	"intake",
	"intake-parquet",
	"ncls",
    # "pypairix",  FIXFIX: this is supplied by the pairix conda package which breaks the conda.recipe
	"pairtools",
	"pandas >=0.25",
	"numpy >=1.16",
	"pyarrow",
	"pysam",
	"tqdm",
]

extras_require = {
	"tests": [
		"pytest"
	],
	"build": [],
	"doc": [
		"sphinx <2.0",
	]
}

extras_require['dev'] = extras_require['tests'] + extras_require['doc']  + ['black', 'isort', 'flake8']


setup(
    name='pore_c',
    version=version,
    license='MPL 2.0',
    description='Tools to process data from an ONT poreC run',
    long_description='',
    author='Matthew Pendleton',
    author_email='Matthew.Pendleton@nanoporetech.com',
    url='https://github.com/nanoporetech/pore_c',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Utilities',
    ],
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
	python_requires=python_requires,
    install_requires=install_requires,
    tests_require=extras_require['tests'],
    extras_require=extras_require,
    entry_points={
        'console_scripts': [
            'pore_c = pore_c.cli:cli',
        ]
    }
)
