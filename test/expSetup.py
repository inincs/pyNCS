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

def build_setup(setupfile = 'test.xml', setups_dir = 'setupfiles/', chips_dir = 'chipfiles/'):
    nsetup = NeuroSetup(
            setups_dir+'test_setuptype.xml',
            setups_dir+setupfile,
            prefix = chips_dir,
            offline=False)
    set_default_biases(nsetup=nsetup, biases_dir='biases/')
    return nsetup

if __name__ == '__main__':
    nsetup = build_setup()
    nsetup.chips['ifslwta'].set_parameter('pinj', 2.82)    

    seq_pop = pyNCS.Population('default', 'Default Population') 
    seq_pop.populate_by_number(nsetup, 'seq', 'excitatory', 60)
    
    seq_pop1= seq_pop[:30]
    seq_pop2= seq_pop[30:]
    
    exc_pop = pyNCS.Population('default', 'Default Population') 
    exc_pop.populate_by_number(nsetup, 'ifslwta', 'excitatory', 60)

    exc_pop1= exc_pop[:30]
    exc_pop2= exc_pop[30:]
               
    mon = nsetup.monitors.import_monitors_otf([exc_pop1, exc_pop2]) 
       
    c1=pyNCS.PConnection(seq_pop1, exc_pop1,'excitatory0', 'random_all2all',{'p':1.0})
    c2=pyNCS.PConnection(seq_pop2, exc_pop2,'excitatory1', 'random_all2all',{'p':.1})
    
    nsetup.chips['ifslwta'].set_parameter('nsynstdw0',.45)
    nsetup.chips['ifslwta'].set_parameter('nsynstdw1',.45)
    stim = seq_pop.soma.spiketrains_inh_poisson(np.array([np.linspace(1,1000,10)]*60), np.linspace(0,1000,10))
    
    out = nsetup.run(stim)
    from pylab import *
    ion()
    pyNCS.monitors.RasterPlot(mon)
    pyNCS.monitors.MeanRatePlot(mon)
