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

# Traits imports
try:
    from enthought.traits.api import *
    from enthought.traits.ui.api import *
except ImportError:
    from traits.api import *
    from traitsui.api import * #This causes pyNCS not to run in headless mode (without C)


class Parameter(HasTraits):
    SignalName = Str('Parameter name')  # Parameter name
    _onlyGui = False  # Flag to only update GUI
    v = Property(Range(-1., 3.3))
    _v = None

    def _get_v(self):
        if self._v is None:
            self._v = self.getValue()
        return self._v

    def _set_v(self, value):
        if not self._onlyGui:
            self.setValue(value)
        else:
            self._v = value

    view = View(Group(Item('SignalName', style='readonly',
                           show_label=False),
                      Item('v', show_label=False, resizable=True),
                      orientation='horizontal'),
                resizable=True,
               )

    def __init__(self, parameters, configurator):
        '''
        Parameter(parameters, configurator)
        parameters: dictionary of parameters and values
        This object is designed to be used with the configurator to set parameters
        '''
        self.param_data = dict(parameters)
        self.configurator = configurator
        # Initialize variable for GUI
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
        self._v = x  # Update gui
        return x