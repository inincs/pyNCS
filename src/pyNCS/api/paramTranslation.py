#-----------------------------------------------------------------------------
# File Name : sympyParams.py
# Purpose:
#
# Author: Emre Neftci
#
# Creation Date : 
# Last Modified : Tue 15 Mar 2011 06:11:41 PM CET
#
# Copyright : (c) 2011
# Licence : GPLv4
#----------------------------------------------------------------------------- 

#chip=pyNCS.Chip('/home/emre/Neuron/software/pyNCS/src/pyNCS/pads_ifslwta.csv')
import numpy as np
from xml.dom import minidom as md
from urllib2 import urlopen
#chip=pyNCS.Chip('/home/emre/Neuron/software/pyNCS/src/pyNCS/pads_ifslwta.csv')
from numpy import log,exp

def log(x):
    if x<0:
        return -1e32
    else:
        return np.log(x)

class params():

    xml_version=.1


    def __init__(self,configurator, xml_filename):
        from sympy import Symbol,exp,log
        self.conf=configurator
        self.cur=dict()

        for s in self.conf.get_param_names():            
            globals()[s]=Symbol(s)
        
        
        self.loadXML(xml_filename)    

        #Some default parameters
        for i,k in self.initial_parameters.items():
            self.cur.update({i:k})

        self.cur.update(self.conf.get_parameters())                     

    def loadXML(self,filename):
        if '://' in filename:
            doc= md.parse( urlopen(filename) )
        else:
            doc = md.parse( filename )
        
        nsetup= doc.childNodes[0]
        if nsetup.tagName!='paramtrans' or float(nsetup.getAttribute('version'))<self.__class__.xml_version:
            raise Exception( 'wrong file-format' )
 
        for n in nsetup.getElementsByTagName('field'):
            if n.getAttribute('id') == 'constants': self._set_constants(n)
            elif n.getAttribute('id') == 'calibration': self._set_calibration(n)
            elif n.getAttribute('id') == 'defaultparameters': self._set_defaultparameters(n)
            elif n.getAttribute('id') == 'translationnames': self._set_translationnames(n)
            elif n.getAttribute('id') == 'translationfunctions': self._set_translationfunctions(n)
            else: raise RuntimeError("Unknown field in provided xml file: "+n.getAttribute('id'))
 
    def _set_constants(self,n):
        from sympy import Symbol,exp,log
        self.constant_parameters=dict() 
        for p in n.getElementsByTagName('parameter'):
            #Considfer defining a function (for esthetic reasons)
            self.constant_parameters[str(p.getAttribute( 'id' ))]=float(eval(p.getAttribute( 'value' )))
        for i in self.constant_parameters.keys():
            globals()[i]=Symbol(i)
            self.cur[i]=self.constant_parameters[i]

    def _set_calibration(self,n):
        from sympy import Symbol,exp,log
        self.calib_parameters=dict() 
        for p in n.getElementsByTagName('parameter'):
            #Considfer defining a function (for esthetic reasons)
            self.calib_parameters[str(p.getAttribute( 'id' ))]=float(eval(p.getAttribute( 'value' )))
        for i in self.calib_parameters.keys():
            globals()[i]=Symbol(i)
            self.cur[i]=self.calib_parameters[i]

    def _set_defaultparameters(self,n):
        from sympy import Symbol,exp,log
        self.initial_parameters=dict() 
        for p in n.getElementsByTagName('parameter'):
            #Considfer defining a function (for esthetic reasons)
            self.initial_parameters[str(p.getAttribute( 'id' ))]=float(eval(p.getAttribute( 'value' )))
        for i in self.initial_parameters.keys():
            globals()[i]=Symbol(i)
            self.cur[i]=self.initial_parameters[i]

    def _set_translationnames(self,n):
        from sympy import Symbol,exp,log
        self.translation_names=dict() 
        for p in n.getElementsByTagName('name'):
            #Considfer defining a function (for esthetic reasons)
            self.translation_names[str(p.getAttribute( 'standard' ))]=p.getAttribute( 'native' )
        for i in self.translation_names.keys():
            globals()[i]=Symbol(i)
            self.cur[i]=self.translation_names[i]

    def _set_translationfunctions(self,n):
        from sympy import Symbol,solve,exp,log
        self.translation_reverse=dict() 
        for p in n.getElementsByTagName('function'):
            #Considfer defining a function (for esthetic reasons)
            self.translation_reverse[str(p.getAttribute( 'standardname' ))]=eval(p.getAttribute( 'body' ))
            
        for i in self.translation_reverse.keys():            
            self.cur[i]=self.translation_reverse[i]

    def translate(self,varname,val):
        if val<1e-32: val=1e-32
        
        self.cur.update({varname:val})
        if varname in self.translation_names.values():
        #Is transformed name: do forward first
            self.update()
    
    def update(self):
        for i in self.translation_names:
            self.cur.update({i:eval(str(self.translation_reverse[i]),globals(),self.cur)})
 
        
def loadBiases(filename):
    h=file(filename,'r')
    d={}
    while h:
        s=h.readline()
        if len(s)==0:
            break
        else:
            n,v=s.strip().split('\t')
            d[n] = float(v)
    h.close()
    return d

if __name__=='__main__':
    import pyNCS
    configurator = pyNCS.ConfAPI.Configurator()
    configurator._readCSV('chipfiles/ifslwta.csv')
    configurator.set_parameters(loadBiases('biases/defaultBiases_ifslwta'))
    p=params(configurator, 'chipfiles/ifslwta_paramtrans.xml')
    p.update()
    