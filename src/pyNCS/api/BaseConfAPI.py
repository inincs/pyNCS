#-----------------------------------------------------------------------------
# Purpose:
#
# Author: <authors name>
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
#ConfAPI
#Biases and mapper
#Api for modules having pyAMDA-like functionality
#Api for modules having pyAEX-like functionality

from types import GeneratorType  # For checking generator type
from contextlib import contextmanager
from lxml import etree
import warnings

class Parameter:
    def __init__(self, parameters, configurator):
        '''
        Parameter(parameters, configurator)
        parameters: dictionary of parameters and values
        This object is designed to be used with the configurator to set parameters
        '''
        self.param_data = dict(parameters)
        self.configurator = configurator
        self.SignalName = self.param_data['SignalName']

    def __str__(self):
        return str(self.param_data)

    def __getNHML__(self):
        '''
        Returns lxml.etree.Element representatoin of this parameter
        '''
        doc = etree.Element('parameter')
        for n, v in self.param_data.items():
            doc.attrib[n] = str(v)
        return doc

    def __parseNHML__(self, doc):
        '''
        Parse xml file or element tree to generate the object
        '''
        if isinstance(doc, str):
            # parse the file
            doc = etree.parse(doc).getroot()
        else:
            # assuming doc is an lxml Element object.
            assert doc.tag == 'parameter'
        self.param_data = dict(doc.attrib)
        for k, v in self.param_data.items():
            try:
                v = float(v)
            except:
                pass
            self.param_data[k] = v
        self.SignalName = self.param_data['SignalName']

    def getValue(self):
        return self.configurator.get_parameter(self.param_data['SignalName'])

    def setValue(self, value):
        x = self.configurator.set_parameter(self.param_data['SignalName'], value)
        return x
    
    v = property(getValue, setValue)

