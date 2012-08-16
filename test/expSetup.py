import pyNCS, pyAex
import pyNCS.pyST as pyST
from pyAex import aextcpclientcom, aexionetcom, mapconf #Aex xio module, Aex client using server, mapper module
import pyAMDA.api

import time,sys,random
import numpy as np
import pylab
from pyAexServer import ServerStarter
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
