# -*- coding: utf-8 -*-
from distutils.core import setup

setup(name='pyNCS',
	version='10.09beta',
	description='Package for NCS experiment',
	author='Sadique Sheik, Emre Neftci, Fabio Stefanini, Giacomo Indiveri',
	author_email='sadique@ini.phys.ethz.ch, emre@ini.phys.ethz.ch, giacomo@ini.phys.ethz.ch, fabios@ini.phys.ethz.ch',
	url='https://github.com/inincs/pyNCS',
	packages = ['pyNCS', 'pyNCS.AerViewer', 'pyNCS.pyST'],
	package_dir={'pyNCS' : 'src/pyNCS'},
     )
