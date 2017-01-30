# -*- coding: utf-8 -*-
from distutils.core import setup
import os

setup(name='pythonNCS',
	version='20170129',
	description='Python Neurormophic Chips and Systems',
	author='Sadique Sheik, Emre Neftci, Fabio Stefanini, Giacomo Indiveri',
	author_email='sadique@ini.phys.ethz.ch, eneftci@uci.edu, giacomo@ini.phys.ethz.ch, fabios@ini.phys.ethz.ch',
	url='https://github.com/inincs/pyNCS',
        download_url='https://github.com/inincs/pyNCS/tarball/stable_20170129',
	packages = ['pyNCS', 'pyNCS.pyST', 'pyNCS.api'],
        package_dir={'pyNCS' : 'src/pyNCS'},
	package_data={'pyNCS' : ['data/*.dtd',
                             'data/chipfiles/*.csv',
                             'data/chipfiles/*.nhml',
                             'data/chipfiles/*.xml']},
     )
