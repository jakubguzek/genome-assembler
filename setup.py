#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import io
import os.path
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext
from typing import Any

from setuptools import find_packages
from setuptools import setup


def read(*names: Any, **kwargs: dict[str, str]):
    with io.open(join(dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")) as fh:  # type: ignore
        return fh.read()


setup(
    name="Genome Assembler",
    version="0.1.2",
    license="MIT",
    description="A CLI genome assembler based on de Bruijn graph method written in Python ",
    long_description="{}\n{}".format(
        re.compile("^.. start-badges.*^.. end-badges", re.M | re.S).sub(
            "", read("README.md")
        ),
        re.sub(":[a-z]+:`~?(.*?)`", r"``\1``", read("CHANGELOG.md")),
    ),
    long_description_content_type="text/markdown",
    author="Jakub J. Guzek",
    author_email="jakub.j.guzek@gmail.com",
    url="file://" + os.path.abspath(dirname(__file__)),
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        "Topic :: Utilities",
        "Private :: Do Not Upload",
    ],
    keywords=["genome" "assembly", "de Bruijn" "python", "graph"],
    python_requires=">=3.10",
    install_requires=[
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
        "biopython==1.80",
        "contourpy==1.0.7",
        "cycler==0.11.0",
        "fonttools==4.38.0",
        "graphviz==0.20.1",
        "kiwisolver==1.4.4",
        "matplotlib==3.6.3",
        "numpy==1.24.1",
        "packaging==23.0",
        "pandas==1.5.3",
        "Pillow==9.4.0",
        "pyparsing==3.0.9",
        "python-dateutil==2.8.2",
        "pytz==2022.7.1",
        "six==1.16.0",
    ],
)
