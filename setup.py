# -*- coding: utf-8 -*-
u"""rsbeams setup script

:copyright: Copyright (c) 2016 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from pykern.pksetup import setup

setup(
    name='rsbeams',
    author='RadiaSoft LLC',
    author_email='pip@radiasoft.net',
    install_requires=[
        'h5py',
        'pykern',
        'pathos',
        'numpy>=1.19.1',
        'scipy',
        'sympy>=1.2',
        'ruamel.yaml'
    ],
    entry_points={
        'console_scripts': ['kinematic=rsbeams.rsstats.kinematic:main'],
    },
    description='Code-agnostic Python utilities for particle beam simulations',
    license='http://www.apache.org/licenses/LICENSE-2.0.html',
    url='https://github.com/radiasoft/rsbeams',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python',
        'Topic :: Utilities',
    ],
)
