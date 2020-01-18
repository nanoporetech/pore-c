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
    version="0.1",
    description="Tools to process data from an ONT poreC run",
    python_requires="==3.7.*,>=3.7.0",
    author="Matthew Pendleton",
    author_email="Matthew.Pendleton@nanoporetech.com",
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
    package_data={},
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
        "pandas==0.25.*,>=0.25.0",
        "pyarrow",
        "pydantic==1.*,>=1.3.0",
        "pyranges==0.*,>=0.0.71",
        "pysam",
        "streamz==0.*,>=0.5.2",
        "tqdm",
    ],
)
