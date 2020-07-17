# -*- coding: utf-8 -*-

# DO NOT EDIT THIS FILE!
# This file has been autogenerated by dephell <3
# https://github.com/dephell/dephell

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = ""

setup(
    long_description=readme,
    name="pore_c",
    version="0.3.0-rc.1",
    description="Tools to process data from an ONT poreC run",
    python_requires="==3.7.*,>=3.7.0",
    author="Eoghan Harrington",
    author_email="Eoghan.Harrington@nanoporetech.com",
    maintainer="Matthew Pendleton",
    maintainer_email="Matthew.Pendleton@nanoporetech.com",
    license="MPL 2.0",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Utilities",
    ],
    entry_points={"console_scripts": ["pore_c = pore_c.cli:cli"]},
    packages=["pore_c"],
    package_dir={"": "."},
    package_data={"pore_c": ["*.profraw", "analyses/*.profraw"]},
    install_requires=[
        "biopython",
        "click==7.0.*,>=7.0.0",
        "cooler==0.8.*,>=0.8.5",
        "cython>=0.29",
        "dask==2.*,>=2.0.0",
        "distributed==2.*,>=2.9.3",
        "intake",
        "intake-parquet",
        "ncls",
        "networkx==2.*,>=2.4.0",
        "numpy>=1.16.1",
        "pairtools",
        "pandas>=1.0.5",
        "pyarrow==1.*,>=1.0.0",
        "pydantic==1.6.1",
        "pyranges==0.0.71",
        "pysam",
        "streamz==0.*,>=0.5.2",
        "tqdm",
    ],
    extras_require={
        "dev": [
            "black==19.*,>=19.3.0",
            "darglint==1.*,>=1.1.0",
            "dephell",
            "flake8==3.*,>=3.7.0",
            "flake8-docstrings==1.*,>=1.3.0",
            "isort[pyproject]==4.*,>=4.3.0",
            "pre-commit==1.*,>=1.17.0",
            "pydocstyle<4.0",
            "pytest",
            "pytype==2020.*,>=2020.4.1",
            "seed-isort-config==1.*,>=1.9.0",
            "sphinx==2.*,>=2.1.0",
            "sphinx-autoapi==1.*,>=1.2.1",
            "sphinx-autodoc-typehints==1.*,>=1.10.3",
            "sphinx-click==2.*,>=2.3.1",
            "sphinx-jsonschema==1.*,>=1.13.0",
            "sphinx-rtd-theme==0.*,>=0.4.3",
            "sphinxcontrib-apidoc==0.*,>=0.3.0",
            "toml==0.*,>=0.10.0",
            "vulture==1.*,>=1.4.0",
        ],
        "doc": ["sphinx==2.*,>=2.1.0"],
        "tests": ["pytest"],
    },
)
