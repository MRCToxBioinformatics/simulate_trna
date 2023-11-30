#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
from __future__ import absolute_import, print_function

# import re
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages, setup


def read(*names, **kwargs):
    """Read description files."""
    path = join(dirname(__file__), *names)
    with open(path, encoding=kwargs.get('encoding', 'utf8')) as fh:
        return fh.read()



long_description = '{}'.format(
    read('README.rst')
    )

setup(
    name='simulatetrna',
    version='0.1.0',
    description='Generic simulation of tRNA-Seq reads',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    license='The Unlicense',
    author='Tom Smith',
    author_email='tss38@cam.ac.uk   ',
    url='https://github.com/MRCToxBioinformatics/simulate_trna',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(i))[0] for i in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: The Unlicense',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Microsoft',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    project_urls={
        'webpage': 'https://github.com/MRCToxBioinformatics/simulate_trna',
        'Issue Tracker': 'https://github.com/MRCToxBioinformatics/simulate_trna/issues',
        },
    keywords=[
        'tRNA', 'simulate', 'RNASeq'
        ],
    python_requires='>=3.7, <4',
    entry_points={
        'console_scripts': [
            'simulatetrnacli1= simulatetrnaproject.cli_int1:main',
            ]
        }
    )
