import brian_no_units
from brian.units import *
from brian.group import *
import time
import scipy.signal as sg
import scipy.stats as stats
import numpy as np
from brian.stdunits import *
from brian.network import run, Network, network_operation
from brian.connections.connection import Connection
from brian.directcontrol import SpikeGeneratorGroup
from brian.clock import defaultclock
from brian.equations import Equations
from brian.neurongroup import NeuronGroup
from brian.monitor import SpikeMonitor, StateMonitor

synapse_trans_dict = {0 : 'Ii',
                      1 : 'Ii',
                      2 : 'Ia',
                      3 : 'Ia'}
Tsim=5*second
defaultclock.dt = 0.0001*second
#Network Parameters:
N_EP1=124;
N_IP1=4;
N = N_EP1+N_IP1
#VLSI Neuron constants
C         = 1e-12*farad    
KappaN    = 0.75;
KappaP    = 0.66;
I0N       = 5e-14*amp;
I0P       = 5e-16*amp;
Ut        = 0.025*volt;
Delta_t   = (Ut/KappaP)

#VLSI neuron parameters
V_t    = 0.75*volt
I_in_e = 0
I_in_i = 0
I_lk_e = 2e-12
I_lk_i = 1e-12

#Synaptic variables
tau        =dict();
tau['e'] = 0.1*second
tau['i'] = 0.05*second        #Inhibitory synapse weight
tau['a'] = 0.05*second        #AER synapse time constant

#Feed-Forward parameters
w        =dict()
w['e']    =0e-11
w['i']    =0e-11/N_IP1
w['ei']    =0e-11/N_EP1
w['a']    =0.22e-11
w['ai']    =0.23e-11
    
def create_netobjs(stim):
    ##Feed-Back parameters
    #w        =dict()
    #w['e']    =1e-11
    #w['i']    =5e-11/N_IP1
    #w['ei']    =2e-11/N_EP1
        
    eqs    =Equations('''
    dV/dt=(-I_lk + Ii)/C+(I_fb+I_in+Ia+Ie)/C: volt
    
    
    I_fb = I0P*exp((V-V_t)/Delta_t) : amp
    I_in  : amp
    I_lk  : amp
    
    dIa/dt=-Ia/tau['a']: amp
    dIe/dt=-Ie/tau['e']: amp
    dIi/dt=-Ii/tau['i']: amp
    ''')
        
    EIP    = NeuronGroup(N, model=eqs, reset = 0, threshold=1.5, refractory=.001)
    EP1    = EIP[:N_EP1]
    IP1    = EIP[N_EP1:]
    
    EP1.I_in = np.zeros(N_EP1)*I_in_e
    IP1.I_in = np.zeros(N_IP1)*I_in_i
    EP1.I_lk = np.zeros(N_EP1)*I_lk_e
    IP1.I_lk = np.zeros(N_IP1)*I_lk_i
    

    
    #Create connections between population
    W=(np.tril(np.ones([N_EP1,N_EP1]),1)*np.triu(np.ones([N_EP1,N_EP1]),-1)-np.eye(N_EP1))*w['e']*amp
    ConnE    =Connection(EP1,EP1,'Ie');ConnE.connect(EP1,EP1,W)
    ConnEI    =Connection(EP1,IP1,'Ie');ConnEI.connect(EP1,IP1,w['ei'])
    ConnIE    =Connection(IP1,EP1,'Ii');ConnIE.connect(IP1,EP1,-w['i'])
    
    M_EIP =SpikeMonitor(EIP)    
    MV_EIP= StateMonitor(EIP,'V',record=range(0,N),timestep=int(1*ms/defaultclock.dt))

    @network_operation
    def update_mpot():
        EIP.V[EIP.V<0.]=0.
#     ME_EP1= StateMonitor(EP1,'Ie',record=range(0,N_EP1),timestep=int(1*ms/defaultclock.dt))
#     MI_EP1= StateMonitor(EP1,'Ii',record=range(0,N_EP1),timestep=int(1*ms/defaultclock.dt))
#     MW_EP1= StateMonitor(EP1,'Ia',record=range(0,N_EP1),timestep=int(1*ms/defaultclock.dt))
    
    netobjs =     {'EIP':EIP,                
                   'update_mpot':update_mpot,                   
                   'connE': ConnE,
                   'ConnEI': ConnEI,
                   'ConnIE': ConnIE,
                   #'M_In1': M_In1,
                   'M_EIP': M_EIP,
                   'MV_EIP': MV_EIP}
                   
    return netobjs, M_EIP, MV_EIP
    


if __name__ == '__main__':
    from pylab import *
    from brian.plotting import raster_plot
    ion()
    
    stim = np.transpose([np.random.randint(0,128,32000), np.cumsum(np.random.random(32000)/16)])
    netobjs,  M_In1, M_EIP, MV_EIP  = create_netobjs(stim)
    net = Network(netobjs)
    net.reinit()
    net.run(1)
    raster_plot(*[M_In1, M_EIP])
    

