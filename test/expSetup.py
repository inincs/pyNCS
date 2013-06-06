import pyNCS.pyST as pyST
import time,sys,random
import pyAex, pyNCS
import numpy as np
import pylab
from pyNCS.neurosetup import NeuroSetup


# C O N F I G # # # # # # # # # # # # # # # # # # # # # #

et=pyNCS.et

# set dirnames

def build_setup(setups_dir = 'setupfiles/', chips_dir = 'chipfiles/'):
    nsetup = NeuroSetup(
            setups_dir+'test_setuptype.xml',
            setups_dir+'test.xml',
            prefix = chips_dir,
            offline=False)
    return nsetup
