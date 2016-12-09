#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Sadique Sheik
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
from __future__ import with_statement, absolute_import
from xml.dom import minidom as md
from . import pyST
from .pyST.STas import addrSpec
from .pyST.STas import _buildGrid
import warnings
from .api.ConfAPI import Configurator
from contextlib import contextmanager
from lxml import etree
from re import search as res


class Block():
    def __init__(self, neurochip):
        '''
        A group of items which share some common properties.
        Can represent soma or synapses.
        Direct initialization not intended.
        Should be obtained from the setup or chip NeuronChip object.
        '''
        self.neurochip = neurochip

    def __generateFromXml__(self, chip, xml):
        '''
        Generates the parameter translation xml file.
        NOTE: Backward compatibility function
        '''
        self.domdoc = xml
        self.id = self.domdoc.getAttribute('id')
        self.parameters = {}
        for i in self.domdoc.getElementsByTagName('parameter'):
            param_id = i.getAttribute('id')
            paramname = i.getAttribute('biasname')
            param = self.neurochip.configurator.parameters[paramname]
            self.parameters[param_id] = param
        #Determine if its soma or synapse
        if self.domdoc.tagName == 'soma':
            self.TYPE = 'SOMA'
        elif self.domdoc.tagName == 'synapse':
            self.TYPE = self.domdoc.getAttribute('type').upper() + ' SYNAPSE'
        # TODO: store addresses list instead of ranges, i.e. [[0, 0], [0, 1],
        # ..., [5, 3], [5, 4], ... ]
        # Determining the dimensions of the Group
        # TODO: incapsulate dimensions in neuron group
        self.dims = {}
        for i in self.domdoc.getElementsByTagName('dim'):
            d = i.getAttribute(
                'id')  # ID here is the coordinate (like X , S etc)
            if not i.getAttribute('populate'):
                if d in self.dims:
                    self.dims[d] = self.dims[d] + eval(i.getAttribute('range'))
                else:
                    self.dims[d] = eval(i.getAttribute('range'))
        return self

    def __getNHML__(self):
        '''
        Returns nhml representation of the object
        '''
        if self.TYPE == 'SOMA':
            doc = etree.Element('soma')
        else:
            doc = etree.Element('synapse')
        doc.attrib['type'] = self.TYPE
        for k, v in self.dims.items():
            dimdoc = etree.Element('dim')
            dimdoc.attrib['id'] = k
            dimdoc.attrib['range'] = str(v)
            doc.append(dimdoc)
        for k, v in self.parameters.items():
            paramdoc = etree.Element('parameter')
            paramdoc.attrib['id'] = k
            paramdoc.attrib['SignalName'] = v.SignalName
            doc.append(paramdoc)
        return doc

    def __parseNHML__(self, doc):
        '''
        Parse nhml representation of the object
        '''
        if isinstance(doc, str):
            # parse the file
            doc = etree.parse(doc).getroot()
        else:
            # assuming doc is an lxml Element object
            assert doc.tag == 'soma' or doc.tag == 'synapse'
        self.TYPE = doc.get('type')
        self.dims = {}
        self.parameters = {}
        self.id = doc.get('id')
        for elm in doc:
            if elm.tag == 'dim':
                if self.dims.has_key(elm.get('id')):
                    self.dims[elm.get('id')] = \
                        self.dims[elm.get('id')] + eval(elm.get('range'))
                else:
                    self.dims[elm.get('id')] = eval(elm.get('range'))

            elif elm.tag == 'parameter':
                self.parameters[elm.get('id')] = \
                        self.neurochip.configurator.parameters[
                            elm.get('SignalName')]
            else:
                pass
        return


class NeuronBlock():
    def __init__(self, neurochip):
        '''
        A specialized 'Group' to encapsulate Neuron populations.
        TODO: No idea what it should be 'specialized' for yet.
        '''
        self.neurochip = neurochip
        self.synapses = {}

    def __generateFromXml__(self, xml):
        '''
        Parse xml file for object representation
        NOTE: Backward compatibility function
        '''
        self.domdoc = xml
        if self.domdoc.tagName == 'neuron':
            self.id = self.domdoc.getAttribute('id')
            # If it is a neuron look for soma and synapses
            self.TYPE = 'NEURON'
            i = self.domdoc.getElementsByTagName(
                'soma')[0]  # there is only one soma per neuron block
            group_id = i.getAttribute('id')
            group = Block(self.neurochip)
            group.__generateFromXml__(self.neurochip, i)
            self.soma = group
            self.soma.addresses = self.expand_dims()

            self.synapses = {}
            for i in self.domdoc.getElementsByTagName('synapse'):
                group_id = i.getAttribute('id')
                group = Block(self.neurochip)
                group.__generateFromXml__(self.neurochip, i)
                self.synapses[group_id] = group
                self.synapses[group_id].addresses = self.expand_dims(group_id)
        else:
            print("This is not a neuron!")

    def __getNHML__(self):
        '''
        Returns nhml representation for this object.
        '''
        doc = etree.Element('neuron')
        # Soma
        somadoc = self.soma.__getNHML__()
        somadoc.attrib['id'] = self.soma.id
        doc.append(somadoc)
        # Synapses
        for k, v in self.synapses.items():
            syndoc = v.__getNHML__()
            syndoc.attrib['id'] = k
            doc.append(syndoc)
        return doc

    def __parseNHML__(self, doc):
        '''
        Parse NHML representation of the object
        '''
        if isinstance(doc, str):
            # parse the file
            doc = etree.parse(doc).getroot()
        else:
            # assuming doc is an lxml Element object
            assert doc.tag == 'neuron'
        self.id = doc.get('id')
        self.synapses = {}
        for elm in doc:
            if elm.tag == 'soma':
                group = Block(self.neurochip)
                group.__parseNHML__(elm)
                self.soma = group
            elif elm.tag == 'synapse':
                group = Block(self.neurochip)
                group.__parseNHML__(elm)
                self.synapses[elm.get('id')] = group
        for k, v in self.synapses.items():
            self.synapses[k].addresses = self.expand_dims(k)
        self.soma.addresses = self.expand_dims()

    def expand_dims(self, synapse=None):
        """
        Returns the list of all addresses for the output (soma) type or the specified input (synapse).
        WARNING: This is buggy! The order of the numbers matters! It should
                 be coherent with the address specification!
        """
        # SADIQUE: Shouldn't this function be modified to use the addrGroups
        # used in populations to be consistent ?
        dims_array = []
        if synapse is None:
            keys = [i['id'] for i in self.neurochip.aerOut.addrConf]
        else:
            keys = [i['id'] for i in self.neurochip.aerIn.addrConf]
        # here they should be in the same order of the
        # address specification
        for k in [i for i in keys if i in self.soma.dims.keys()]:
            dims_array.append(self.soma.dims[k])
        if synapse:
            s = self.synapses[synapse]
            for k in [i for i in keys if i in s.dims.keys()]:
                dims_array = dims_array
                dims_array.append(s.dims[k])  # do it with dict.values()!
        if not len(dims_array) > 0:
            raise Exception('Block has no dimension.')
        return _buildGrid(dims_array)


class Chip:
    def __init__(self, chipdoc, id='noname', offline=False,
                 conf_api=None, conf_kwargs={}):
        '''
        Chip class emulates a chip in the neuromorphic setup
        chipdoc: csv or xml description of the chip
        id      : name of the chip
        offline : True|Fase operation mode.
        conf_api: API library to used for initializing the configurator
        conf_kwargs: dict with all the arguments to the passed on to the
        configurator constructor
        '''
        self.id = id
        # By default the chip is virtual since we do not know the configurator
        # yet.
        self._virtual = True
        self.configurator = Configurator()

        self.chipclass = ''
        # Initialize configurator
        if not offline and conf_api is not None:
            try:
                self.configurator = conf_api.Configurator(**conf_kwargs)
                self.virtual = False
            except Exception as e:
                self.virtual = True
                #Chip should be non-virtual, warn user
                warnings.warn(
                    'Failed to initialize configurator, building offline')
                warnings.warn('Error is :' + str(e))
        else:
            #Chip is intended to be virtual, no need to warn user
            self.virtual = True
        # Parse chip file
        self.dims = []
        if isinstance(chipdoc, str):
            # File name
            if chipdoc.endswith('.csv'):
                self._readCSV(chipdoc)
            elif chipdoc.endswith('.nhml'):
                self.__parseNHML__(chipdoc)
        elif isinstance(chipdoc, etree._Element):
            assert chipdoc.tag == 'chip'
            self.__parseNHML__(chipdoc)
        else:
            pass

    @property
    def virtual(self):
        return self._virtual

    @virtual.setter
    def virtual(self, value):
        if value == False:
            self._virtual = value
        elif value == True:
            self.configurator = Configurator()
            self._virtual = value
        else:
            raise ValueError('virtual value must be a boolean')

    def _readCSV(self, CSVfile):
        '''
        Parse the CSV file to build the chip object.
        NOTE: Backward compatibility function
        '''
        #The following builds the input output address specification only
        self.aerIn, self.aerOut = pyST.STas.load_stas_from_csv(CSVfile)
        try:
            for d in self.aerOut.addrConf:
                self.dims.append(len(d['range']))
        except AttributeError as e:
            warnings.warn('Dimensions could not be inferred')
            warnings.warn('Error is :' + str(e))
        ##The following finds chipclass only
        with open(CSVfile, 'r') as CSV:
            csv = CSV.readlines()
        for line in csv:
            line = line.replace('\'', '')  # remove single and double quotes
            line = line.replace('"', '')
            if 'chipclass' in line.lower():  # Chipclass
                self.chipclass = line.strip().split('\t')[1]
                break
        # Update configurator
        self.configurator._readCSV(CSVfile)
        return

    def __getNHML__(self):
        '''
        Returns nhml representatoin of this object
        '''
        doc = etree.Element('chip')
        doc.attrib['chipclass'] = self.chipclass
        try:
            # aerIn
            aerInDoc = self.aerIn.__getNHML__()
            aerInDoc.attrib['type'] = 'aerIn'
            doc.append(aerInDoc)
        except AttributeError as e:
            warnings.warn("Cannot retreive aerIn data")
        try:
            # aerOut
            aerOutDoc = self.aerOut.__getNHML__()
            aerOutDoc.attrib['type'] = 'aerOut'
            doc.append(aerOutDoc)
        except AttributeError as e:
            warnings.warn("Cannot retreive aerOut in data")
        # Parameters
        doc.append(self.configurator.__getNHML__())
        return doc

    def __parseNHML__(self, doc):
        '''
        Parse nhml file or lxml Element tree to build the object
        '''
        if isinstance(doc, str):
            # parse the file
            doc = etree.parse(doc).getroot()
        else:
            # assuming doc is an lxml Element object
            assert doc.tag == 'chip'
        # Chipclass
        self.chipclass = doc.get('chipclass')
        self.aerIn, self.aerOut = pyST.STas.load_stas_from_nhml(doc)
        for elm in doc:
            if elm.tag == 'parameters':
                self.configurator.__parseNHML__(elm)
            else:
                pass
        # Infer dimensions of the chip
        try:
            for d in self.aerOut.addrConf:
                self.dims.append(len(d['range']))
        except AttributeError as e:
            warnings.warn('Dimensions could not be inferred')
            warnings.warn('Error is :' + str(e))

    @contextmanager
    def context_open(self):
        '''
        Context for opening configurator only if necessary, and closing it if it had to be opened - leaves configurator in its previous state
        '''
        close_when_done = False
        if not self.configurator.isopen:
            self.configurator.open()
            close_when_done = True
        yield
        if close_when_done:
            self.configurator.close()

    def restart(self):
        '''
        Restart the chip
        '''
        with self.context_open():
            self.configurator.reset()

    def get_parameters(self, param_names=None):
        '''
        Get bias values
        '''
        with self.context_open():
            r = self.configurator.get_parameters(param_names)
        return r

    def set_parameters(self, param_dict):
        '''
        Set biases
        '''
        with self.context_open():
            r = self.configurator.set_parameters(param_dict)
        return r

    def set_parameter(self, param_name, param_value):
        '''
        Set biases
        '''
        with self.context_open():
            r = self.configurator.set_parameter(param_name, param_value)
        return r

    def get_parameter(self, param_name):
        '''
        Set biases
        '''
        with self.context_open():
            r = self.configurator.get_parameter(param_name)
        return r

    def get_param_names(self):
        '''
        Returns a list of all bias names
        '''
        return self.configurator.get_param_names()

    def grep_params(self, string):
        '''
        Lists all parameters and values containing the given string.
        Regular expressions allowed, e.g.:
                chip.grep_params('synaer[i]*|synstd[e]*')
        '''
        l = []
        for n in self.get_param_names():
                if res(string, n):
                        l.append(n)
        return self.get_parameters(l)

    def load_parameters(self, CSVfile, sep='\t'):
        '''
        loadBiases(CSVfile) loads the bias valuess from the CSV file
        created with saveBiases(CSVfile)
        '''
        b = {}
        with open(CSVfile, 'r') as CSV:
            #CSV=file(CSVfile)
            for line in CSV.readlines():
                if line.startswith('#') or line.startswith('\n') or line.startswith('\t'):
                    #for commented lines in bias file
                    pass
                else:
                    try:
                        line = line[:line.index('#')]
                    except:
                        pass
                    x = line.strip().split(sep)
                    b[x[0]] = float(x[1])
        self.set_parameters(b)
        print('Parameters loaded and set')
        return None

    def save_parameters(self, filename, *kwargs):
        '''
        Saves all the biases of the chip to a file
        '''
        self.configurator.save_parameters(filename, *kwargs)


class NeuroChip(Chip):
    def __init__(self, chipfile, id='noname', offline=False,
                 conf_api=None, conf_kwargs={}):
        '''
        NeuronChip(chipfile, id='noname',
            configurator=pyAMDA.api.Configurator,
            conf_kwargs={'host':'localhost','board':'206'})
        '''
        self.neuron = {}
        if chipfile.endswith('.csv'):
            chipfile = chipfile.rpartition('.')[0]
            CSVfile = chipfile + ".csv"
            paramfile = chipfile + ".xml"
            Chip.__init__(self, CSVfile, id=id, offline=offline,
                          conf_api=conf_api,
                          conf_kwargs=conf_kwargs)
            self._load_paramfile(paramfile)
        elif chipfile.endswith('.nhml'):
            Chip.__init__(self, chipfile, id=id, offline=offline,
                          conf_api=conf_api,
                          conf_kwargs=conf_kwargs)
            self.__parseNHML__(chipfile)
        else:
            raise Exception('Unknown file format for chipfile!')

    def _load_paramfile(self, paramfile):
        '''
        Function to parse the xml file with the parameter details of the chip.
        '''
        xml = md.parse(paramfile)
        neuron_v = xml.getElementsByTagName('neuron')
        if len(neuron_v):
            for i in neuron_v:
                block_id = i.getAttribute('id')
                neuronblock = NeuronBlock(self)
                neuronblock.__generateFromXml__(i)
                self.neuron[block_id] = neuronblock
        self.domdoc = xml

    def __getNHML__(self):
        '''
        Returns nhml representation of this object.
        '''
        doc = Chip.__getNHML__(self)
        for k, v in self.neuron.items():
            neurodoc = v.__getNHML__()
            neurodoc.attrib['id'] = k
            doc.append(neurodoc)
        return doc

    def __parseNHML__(self, doc):
        '''
        Parse nhml file or lxml Element tree to build the object
        '''
        Chip.__parseNHML__(self, doc)
        if isinstance(doc, str):
            # parse the file
            doc = etree.parse(doc).getroot()
        else:
            # assuming doc is an lxml Element object
            assert doc.tag == 'chip'
        for elm in doc:
            if elm.tag == 'neuron':
                neuronblk = NeuronBlock(self)
                neuronblk.__parseNHML__(elm)
                self.neuron[elm.get('id')] = neuronblk
            else:
                pass
        # Infer dimensions of the chip
        try:
            for d in self.aerOut.addrConf:
                self.dims.append(len(d['range']))
        except AttributeError as e:
            warnings.warn('Dimensions could not be inferred')
            warnings.warn('Error is :' + str(e))

    def saveNHML(self, filename):
        '''
        save NHML representation of this object to a file
        '''
        doc = self.__getNHML__()
        with open(filename, 'w') as f:
            f.write(etree.tostring(doc, pretty_print=True))
