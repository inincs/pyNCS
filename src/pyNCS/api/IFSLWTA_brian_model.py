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


defaultclock.dt = 0.0001*second
#Network Parameters:
N_EP1=124;
N_IP1=4;
N = N_EP1+N_IP1
#VLSI Neuron constants
sigma_mismatch=0.25
    
def create_netobjs(stim, params):
    C         = params['Cap']   
    KappaN    = params['Kappan']
    KappaP    = params['Kappap']
    I0N       = params['I0n']
    I0P       = params['I0p']
    Ut        = params['Ut']
    Delta_t   = (Ut/KappaP)
    
    #Feed-Forward parameters
    i_inj           =       params['i_inj'    ]
    i_injinh1       =       params['i_injinh1']
    i_injinh2       =       params['i_injinh2']
    i_injinh3       =       params['i_injinh3']
    i_injinh4       =       params['i_injinh4']
    i_leak          =       params['i_leak'   ]
    tau_syn_E       =       params['tau_syn_E']
    tau_syn_I       =       params['tau_syn_I']
    tau_synloc_E    =       params['tau_synloc_E']
    tau_synloc_IE   =       params['tau_synloc_IE']
    v_thresh        =       params['v_thresh' ]
    w_syn_E1        =       params['w_syn_E1' ]
    w_syn_E2        =       params['w_syn_E2' ]
    w_syn_I         =       params['w_syn_I'  ]
    w_synloc_E      =       params['w_synloc_E']
    w_synloc_EE     =       params['w_synloc_EE']
    w_synloc_EI     =       params['w_synloc_EI']
    w_synloc_IE     =       params['w_synloc_IE']
    w_synloc_S      =       params['w_synloc_S']
    ##Feed-Back parameters
    #w        =dict()
    #w['e']    =1e-11
    #w['i']    =5e-11/N_IP1
    #w['ei']    =2e-11/N_EP1
        
    eqs    =Equations('''
    dV/dt=(-I_lk + I_fb + I_in + Ia + Iloce - Iloci - Ii)/C: volt
    
    
    I_fb = I0P*exp((V-v_thresh)/Delta_t) : amp
    I_in : amp
    I_lk = i_leak : amp
    
    dIa/dt=-Ia/tau_syn_E: amp
    dIloce/dt=-Iloce/tau_synloc_E: amp
    dIloci/dt=-Iloci/tau_synloc_IE: amp
    dIi/dt=-Ii/tau_syn_I: amp      
    ''')
        
    EIP    = NeuronGroup(N, model=eqs, reset = 0, threshold=1.5, refractory=.001)
    EP1    = EIP[:N_EP1]
    IP1    = EIP[N_EP1:]
    
    EP1.I_in = np.random.normal(1,sigma_mismatch,N_EP1)*i_inj
    IP1.I_in = np.random.normal(1,sigma_mismatch,N_IP1)*\
                                   np.array([i_injinh1,
                                             i_injinh2,
                                             i_injinh3,
                                             i_injinh4])
    
    #Create connections between population
    v_loclat = np.zeros([N_EP1])
    v_loclat[N_EP1/2] = w_synloc_S
    v_loclat[[N_EP1/2-1,N_EP1/2+1]] = w_synloc_E
    v_loclat[[N_EP1/2-2,N_EP1/2+2]] = w_synloc_EE
    v_loclat[[N_EP1/2-3,N_EP1/2+3]] = w_synloc_EE/2
    v_loclat = np.roll(v_loclat,-N_EP1/2)
    W = np.array([ np.roll(v_loclat,i) for i in range(N_EP1)])
    W *= np.random.normal(1,sigma_mismatch,W.shape)    
    ConnE    =Connection(EP1,EP1,'Iloce');ConnE.connect(EP1,EP1,W)
    ConnEI   =Connection(EP1,IP1,'Iloce');ConnEI.connect(EP1,IP1, W = w_synloc_EI*np.random.normal(1,sigma_mismatch,[len(EP1),len(IP1)]))
    ConnIE   =Connection(IP1,EP1,'Iloci');ConnIE.connect(IP1,EP1, W = w_synloc_IE*np.random.normal(1,sigma_mismatch,[len(IP1),len(EP1)]))
    
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
                   'ConnE': ConnE,
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
    from paramTranslation import params, loadBiases
    from expSetup import *
    configurator = pyNCS.ConfAPI.Configurator()
    configurator._readCSV('chipfiles/ifslwta.csv')
    configurator.set_parameters(loadBiases('biases/defaultBiases_ifslwta'))
    p=params(configurator, 'chipfiles/ifslwta_paramtrans.xml')
    p.translate('pinj',2.8)    
    p.translate('nsynloclat1',0.55)
    p.translate('nsynloclat2',0.53)
    p.translate('nsynexcinh',0.45)   
    p.translate('psynlocinhw',2.75)
    p.translate('psynlocinhth',.3)
    stim = np.transpose([np.random.randint(0,128,32000), np.cumsum(np.random.random(32000)/16)])
    netobjs,  M_EIP, MV_EIP  = create_netobjs(stim,p.cur)
    net = Network(netobjs.values())
    net.run(1)
    raster_plot(*[M_EIP])
    


