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
setupdir = 'setupfiles/'

def build_setup():
    nsetup = NeuroSetup(
            setupdir+'test_setuptype.xml',
            setupdir+'test.xml',
            offline=False)
    return nsetup
