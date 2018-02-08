# -*- coding: utf-8 -*-
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy as np

setup(name='pythonNCS',
	version='20170129',
	description='Python Neurormophic Chips and Systems',
	author='Sadique Sheik, Emre Neftci, Fabio Stefanini, Giacomo Indiveri',
	author_email='sadique@ini.phys.ethz.ch, eneftci@uci.edu, giacomo@ini.phys.ethz.ch, fabios@ini.phys.ethz.ch',
	url='https://github.com/inincs/pyNCS',
        download_url='https://github.com/inincs/pyNCS/tarball/stable_20170129',
	include_dirs=[np.get_include()],
	ext_modules = cythonize('src/pyNCS/pyST/*.pyx'),
	packages = ['pyNCS', 'pyNCS.pyST', 'pyNCS.api'],
        package_dir={'pyNCS' : 'src/pyNCS'},
	package_data={'pyNCS' : ['data/*.dtd',
                             'data/chipfiles/*.csv',
                             'data/chipfiles/*.nhml',
                             'data/chipfiles/*.xml']},
     )
