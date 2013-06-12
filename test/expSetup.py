import pyNCS.pyST as pyST
import time,sys,random
import pyAex, pyNCS
import numpy as np
import pylab
from pyNCS.neurosetup import NeuroSetup
import warnings

# C O N F I G # # # # # # # # # # # # # # # # # # # # # #

et=pyNCS.et

# set dirnames
def set_default_biases(nsetup=None, biases_dir='biases/'):
    for c in nsetup.chips.itervalues():

        filename=''
        try:
            if not c.virtual:
                filename=biases_dir+'defaultBiases_'+c.chipclass.lower()
                c.load_parameters(filename)
        except IOError as (errno, strerror):
            warnings.warn("I/O error({0}): {1}".format(errno, strerror))
            warnings.warn("Could not find file {0}".format(filename))
            pass

def build_setup(setups_dir = 'setupfiles/', chips_dir = 'chipfiles/'):
    nsetup = NeuroSetup(
            setups_dir+'test_setuptype.xml',
            setups_dir+'test.xml',
            prefix = chips_dir,
            offline=False)
    set_default_biases(nsetup=nsetup, biases_dir='biases/')
    return nsetup
