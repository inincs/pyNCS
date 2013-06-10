#!/bin/python
#-----------------------------------------------------------------------------
# File Name : testComAPI.py
# Purpose: Brian Simulator API for testing pyNCS. Can only be used as one 
# instance because mappings are global 
#
# Author: Emre Neftci
#
# Creation Date : 05-06-2013
# Last Modified : Wed 05 Jun 2013 06:47:18 PM PDT
#
# Copyright : (c) 
# Licence : GPLv2
#----------------------------------------------------------------------------- 

from pyNCS.ComAPI import *
from pyNCS.ConfAPI import *
import IFSLWTA_brian_model as ibm
import numpy as np



chipname = 'ifslwta'
channel_seq = 0
channel_ifslwtain = 6
global_mappings_list = []
global_neurosetup = None

def _decode_mappings_list(mappings_list, synapse_id):
    ml = np.array(mappings_list)
    fr = decode_addr(ml[:,0])[channel_seq][0] #0 for field 0
    to = decode_addr(ml[:,1])[channel_ifslwtain]
    to_syn = to[0,to[1,:]==synapse_id]
    if len(to_syn)==0:
        return None
    else:        
        if np.shape(ml)[1]==3:
            prob = ml[:,2].astype('float')/128
        else:    
            prob = np.ones(ml.shape[0])
        return zip(fr,to_syn,prob)

def _dlist_to_mappingdict(mapping):
    from collections import defaultdict
    #sort list
    mapping_dict = defaultdict(list)
    mapping_dict_probs = defaultdict(list)
    def func(srctgt):
        mapping_dict[srctgt[0]].append(srctgt[1])
        mapping_dict_probs[srctgt[0]].append(srctgt[2])
    map(func, mapping)
    return mapping_dict, mapping_dict_probs

def _mappingdict_to_matrix(mapping_dict, mapping_dict_probs):
    N_fr = len(mapping_dict)
    N_to = max(max(mapping_dict.values()))+1
    M = np.zeros([N_fr, N_to])
    P = np.zeros([N_fr, N_to])
    for k in mapping_dict.keys():
        M[k,mapping_dict[k]] = 1
        P[k,mapping_dict[k]] = mapping_dict_probs[k] 
    return M, P

def translate_mappings(mapping_list, synapse_id):
    m = _decode_mappings_list(mapping_list, synapse_id)
    if m ==None:
        return None
    else:
        return _mappingdict_to_matrix(*_dlist_to_mappingdict(m))

def decode_addr(addr):
    #decode physical addresses to neuron - synapse    
    global global_neurosetup    
    return global_neurosetup.seq.addrPhysicalExtract(addr)    

class Communicator(BatchCommunicatorBase):
    def run(self, stimulus=None, duration=None, context_manager=None):
        stimulus_abs = stimulus.copy().astype('float')
        stimulus_abs[:,1] = np.cumsum(stimulus_abs[:,1])/1000000
        if duration == None:
            duration = 1
        elif duration == None and stimulus is not None:
            duration = np.max(stimulus_abs[:,1])/1000
        else:
            duration = float(duration)/1000
        from brian.network import Network        
        net, M_EIP, MV_EIP = self._prepare_brian_net(stimulus_abs)
        net.reinit(states=False)   
        net.run(duration)
        sp = np.array(M_EIP.spikes)
        sp[:,0]+=(2<<13) #slot 2 in text.xml TODO
        sp[:,1]*=1000 #slot 2 in text.xml TODO
        return sp
        
    
    def _prepare_brian_net(self, stimulus):
        from brian.directcontrol import SpikeGeneratorGroup
        from brian.synapses.synapses import Synapses
        from brian.network import Network

        netobjs, M_EIP, MV_EIP = ibm.create_netobjs(stimulus)
        gml = np.array(global_mappings_list)
        N = max(gml[:,0])+1
        W,P = translate_mappings(global_mappings_list, 2)
        In1=SpikeGeneratorGroup(N,stimulus[stimulus[:,0] <N])
        stimulus[:,1]=stimulus[:,1]
        stimulus=stimulus[stimulus[:,1].argsort()]
        N = max(stimulus[:,0])    
        EIP = netobjs['EIP']
        S=Synapses(In1,EIP,model="""w : 1
                                    p : 1""",
                         pre="Ia+=w*(rand()<p)")
        gml_dict = _dlist_to_mappingdict(_decode_mappings_list(gml,2)) 
        for i in gml_dict[0]:
            S[i,gml_dict[0][i]]=True
            S.w[i,gml_dict[0][i]] = ibm.w['a']
            S.p[i,gml_dict[0][i]] = gml_dict[1][i]  
        self.S = S
        net = Network(netobjs.values()+[S,In1])
        return net, M_EIP, MV_EIP  
        
    
class Configurator(ConfiguratorBase):
    def register_neurosetup(self, neurosetup):
        '''
        Provides a link to the Neurosetup. This is useful for complex parameter
        configuration protocols requiring the sequencing and monitoring of
        address-events
        '''
        self._neurosetup_registered = True
        self._neurosetup = neurosetup
        global global_neurosetup
        global_neurosetup = neurosetup
        
    def get_parameter(self, param_name):
        #IMPLEMENT
        '''
        Gets parameter param_name.
        '''
        return 0
    
    def set_parameter(self, param_name, param_value):
        #IMPLEMENT
        '''
        Sets parameter param_name with param_value
        '''
        return None
    
    def reset(self):
        #IMPLEMENT
        '''
        Resets all the parameters to default values
        '''
        return None
    
class Mappings(MappingsBase):
    def add_mappings(self, mappings):
        #IMPLEMENT (REQUIRED)
        '''
        Adds *mappings* to the mappings table.

        Inputs:
        *mappings*: a two-dimenstional iterable
        '''
        global global_mappings_list
        global_mappings_list += mappings 

    def get_mappings(self):
        #IMPLEMENT (REQUIRED)
        '''
        Returns an array representing the mappings
        '''
        global global_mappings_list        
        return global_mappings_list

    def clear_mappings(self):
        #IMPLEMENT (REQUIRED)
        '''
        Clears the mapping table. No inputs
        '''
        global global_mappings_list
        while True:
            try:
                global_mappings_list.pop()
            except IndexError:
                break
            
        return None

