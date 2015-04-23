#!/bin/python
#-----------------------------------------------------------------------------
# File Name : testComAPI.py
# Purpose: Brian Simulator API for testing pyNCS. Can only be used as one 
# instance because mappings are global 
#
# Author: Emre Neftci
#
# Creation Date : 05-06-2013
# Last Modified : Thu 26 Mar 2015 12:21:07 PM PDT
#
# Copyright : (c) 
# Licence : GPLv2
#----------------------------------------------------------------------------- 
from pyNCS.neurosetup import get_data
from pyNCS.api.ComAPI import *
from pyNCS.api.ConfAPI import *
import pyNCS.pyST as pyST
import IFSLWTA_brian_model as ibm
import numpy as np
from brian.directcontrol import PoissonGroup
import os

#TODO: multichip mappings
#TODO: generalize _decode_mappings_list

chipname = 'ifslwta'
channel_seq = 0
channel_ifslwtain = 6
global_mappings_list = []
global_neurosetup = None
global_sympy_params = None

synapse_trans_dict = {0 : ['Ii', 'w_syn_I'],
                      1 : ['Ii', 'w_syn_I'],
                      2 : ['Ia', 'w_syn_E1'],
                      3 : ['Ia', 'w_syn_E2']}

def _decode_mappings_list(mappings_list, synapse_id):
    ml = np.array(mappings_list)
    fr = decode_addr(ml[:,0])[channel_seq] #0 for field 0
    to = decode_addr(ml[:,1])[channel_ifslwtain]
    to_syn = to[0,to[1,:]==synapse_id]
    fr_syn = fr[0,to[1,:]==synapse_id]
    if len(to_syn)==0:
        return []
    else:        
        if np.shape(ml)[1]==3:
            prob = ml[:,2].astype('float')/128
            prob_syn = prob[to[1,:]==synapse_id]
        else:    
            prob = np.ones(ml.shape[0])            
        return zip(fr_syn,to_syn,prob_syn)

def _dlist_to_mappingdict(mapping):
    from collections import defaultdict
    #sort list
    if len(mapping)>0:
        mapping_dict = defaultdict(list)
        mapping_dict_probs = defaultdict(list)
        def func(srctgt):
            mapping_dict[srctgt[0]].append(srctgt[1])
            mapping_dict_probs[srctgt[0]].append(srctgt[2])
        map(func, mapping)
        return mapping_dict, mapping_dict_probs
    else:
        return {},{}

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
        stimulus_abs = pyST.events(stimulus, atype='p', isISI=True)
        stimulus_abs.set_abs_tm()
        self.evn_sin = evs_in = stimulus_abs.get_adtmev().astype('float')
        evs_in[:,1]*=1e-6 #brian takes seconds
        if duration == None and len(stimulus_abs)>0:
            duration = np.max(evs_in[:,1])
            print duration
        elif duration == None:
            duration = 1.
        else:
            duration = float(duration)/1e3
        net, M_EIP, MV_EIP = self._prepare_brian_net(evs_in)
        self.outs = [net, M_EIP, MV_EIP]
        net.reinit(states=False)
        print('running virtual IFLSWTA for {0}s'.format(duration))  
        net.run(duration)
        sp = np.array(M_EIP.spikes).reshape(-1,2)
        sp[:,0]+=(2**16) #slot 2 in text.xml. THIS IS A POSSIBLE SOURCE OF TEST ERROR
        sp[:,1]*=1e6 
        return sp
        
    
    def _prepare_brian_net(self, stimulus):
        from brian.directcontrol import SpikeGeneratorGroup
        from brian.synapses.synapses import Synapses
        from brian.network import Network
        self.S = []
        netobjs, M_EIP, MV_EIP = ibm.create_netobjs(stimulus, global_sympy_params.cur)
        
        if len(stimulus)>0: 
            N = max(stimulus[:,0])+1
            In1=SpikeGeneratorGroup(N,stimulus[stimulus[:,0] <N])
            self.S.append(In1)
        #Create connections, if any
        if len(global_mappings_list)>0:                 
            gml = np.array(global_mappings_list).reshape(-1,3)
            for syn_idx in range(4):   
                      
                gml_dict, pgml_dict = _dlist_to_mappingdict(_decode_mappings_list(gml,syn_idx))
                if len(gml_dict)>0:                
                    if len(stimulus)==0:
                        input_pop = PoissonGroup(max(gml_dict.keys())+1, 0)
                    else:
                        input_pop = In1
                    iname, wname = synapse_trans_dict[syn_idx] 
                    EIP = netobjs['EIP']
                    S=Synapses(input_pop,EIP,model="""w : 1
                                                p : 1""",
                                     pre=iname+"+=w*(rand()<p)")
                     
                    for i in gml_dict:
                        S[i,gml_dict[i]]=True                
                        S.w[i,gml_dict[i]] = global_sympy_params.cur[wname]*np.random.normal(1,ibm.sigma_mismatch, len(gml_dict[i]))
                        S.p[i,gml_dict[i]] = pgml_dict[i] 
                                    
                    self.S.append(S)
            
        net = Network(netobjs.values()+self.S)
                        
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
        return self.parameters[param_name]        
    
    def set_parameter(self, param_name, param_value):
        #IMPLEMENT
        '''
        Sets parameter param_name with param_value
        '''
        param_dict = {param_name:param_value}
        self.sympy_params.cur.update(param_dict)   
        self.sympy_params.update()    
        self.parameters.update(param_dict)
        return None
    
    def set_parameters(self, param_dict):
        #IMPLEMENT
        self.sympy_params.cur.update(param_dict)  
        self.sympy_params.update()     
        self.parameters.update(param_dict)
        return None
    
    
    def add_parameter(self, param):
        #CONVENIENCE FUNCITON. IMPLEMENTATION NOT REQUIRED
        '''
        Add a parameter to the configurator
        param: dictionary with all attributes of parameter or an xml element or
        file name with the <parameter /> element
        '''
        if isinstance(param, dict):            
            self.parameters[param['SignalName']] = 0
        elif isinstance(param, etree._Element):
            self.parameters[param.SignalName] = 0
    
    def reset(self):
        #IMPLEMENT
        '''
        Resets all the parameters to default values
        '''
        return None
    
    def __parseNHML__(self, doc):
        '''
        Parse xml file or element tree to generate the object
        '''
        super(Configurator,self).__parseNHML__(doc)
        global global_sympy_params
        from paramTranslation import params
        filename = get_data('chipfiles/ifslwta_paramtrans.xml')
        global_sympy_params = self.sympy_params = params(self, filename)
        
    def _readCSV(self, CSVfile):        
        super(Configurator,self)._readCSV(CSVfile)
        global global_sympy_params
        from paramTranslation import params
        filename = get_data('chipfiles/ifslwta_paramtrans.xml')
        global_sympy_params = self.sympy_params = params(self, filename)
    
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

