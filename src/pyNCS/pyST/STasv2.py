# -*- coding: utf-8 -*-
# Author: Sadique Sheik
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------

from __future__ import absolute_import

import numpy as np
from lxml import etree
from layoutFieldEncoder import LayoutFieldEncoder
from addressEncoder import AddressEncoder, extract_id_list

class AddrSpec:
    """
    Address specification class
    """
    def __init__(self, doc=None, nhml=True):
        '''
        AddrSpec version
        '''
        self.addrConf = []  
        self.addrPinConf = []  
        self.addrStr = ''  
        self.subAddrSpec = []
        if doc != None:
            self.__parseNHML__(doc)


    def __parseNHML__(self, doc):
        '''
        Parses an lxml element tree or a file name with xml content to
        initialize the object.
        TOexample file
        '''
        if isinstance(doc, str):
            # parse the file
            doc = etree.parse(doc).getroot()
        else:
            # assuming doc is an lxml Element object.
            assert doc.tag == 'addressSpecification'
        
        self.type = doc.get('type')
        self.offset = eval(doc.get('offset'))
        self.min = eval(doc.get('min'))
        self.max = eval(doc.get('max'))

        # addrConf
        for elm in doc:
            if elm.tag == 'pinlayout':
                # Pinlayout
                self.addrStr = elm.text
            elif elm.tag == 'pin':
                # Pins
                pin = {}
                # Id
                pin['id'] = elm.get('id')
                for chld in elm:
                    # Decoder
                    if chld.tag == 'decoder':
                        pin['f'] = chld.text
                    else:
                        pass
                self.addrPinConf.append(pin)
            elif elm.tag == 'dim':
                # Dimensions
                dim = {}
                # Id
                dim['id'] = elm.get('id')
                # Type
                if elm.get('type') == 'synapse':
                    dim['type'] = -1
                elif elm.get('type') == 'soma':
                    dim['type'] = 1
                elif elm.get('type') == 'connection':
                    dim['type'] = -2
                else:
                    # If there are any new types they should be added here!
                    dim['type'] = None
                for chld in elm:
                    if chld.tag == 'range':
                        # Range
                        dim['range'] = eval(chld.text)
                    elif chld.tag == 'default':
                        # Default value
                        dim['default_value'] = eval(chld.text)
                    elif chld.tag == 'description':
                        # Description
                        dim['description'] = chld.text
                    elif chld.tag == 'decoder':
                        # Decoder
                        dim['f'] = chld.text
                    else:
                        pass
                self.addrConf.append(dim)
            elif elm.tag == 'addressSpecification':
                addrSpec = AddrSpec(elm)
                self.subAddrSpec.append(addrSpec)
                ids = extract_id_list(self.addrPinConf)
                try:
                    ids.index(elm.get('id'))
                except ValueError as e:
                    # This field has not been populated yet.
                    self.addrPinConf.append({'id' : elm.get('id'),
                                             'f'  : elm.get('id').lower()})
                    self.addrConf.append({'f' : elm.get('id'),
                                          'id' : elm.get('id').lower(),
                                          'description' : 'Lowlevel physical address space',
                                          'range' : addrSpec.offset + np.arange(addrSpec.min,
                                                                                addrSpec.max)})
            else:
                pass
        # Generate address specification properties
        self.update()
        return

    def update(self):
        '''
        (Re)Generates the object based on the addrStr, addrConf and addrPinConf
        '''
        # Process and complete addrConf object
        self.addrConf = self._process_addrConf()
        
        # Process the no. of dimensions of this address specification object
        self.nDims = len(self.addrConf)
        
        # Generate dictionary holding the values of index of every dimension in
        # addrConf
        self.addrDict = self._process_addrDict()
        
        # Process the bits from the addrStr
        self.addrSpec, self.nBits, self.nBitsTotal = self._stas_parse_addrstr()
        
        # Generate the extract and create functions for physical addresses
        self.field, self.nFields = self._stas_create_fields()
        
        self.addr_encoder = AddressEncoder(self.addrConf, 
                                           self.addrSpec,
                                           self.nBits,
                                           self.nBitsTotal,
                                           self.addrPinConf)
        
        #self.nbits = _stas_compute_nbits(self.addrConf)
        
        #self.addrExtractLogicalFast = dict()
        #self.addrExtractPhysicalFast = dict()


    def _process_addrDict(self):
        '''
        Generate a dictionary holding the values of index of every dimension in addrConf
        '''
        self.addrDict = dict()
        # Find the index of a dimension in addrConf
        for i, hrf in enumerate(self.addrConf):
            self.addrDict[hrf['id']] = i
        return self.addrDict

    def _process_addrConf(self):
        '''
        Find number of required bits from range data
        '''
        for hrf in self.addrConf:
            bits_req = 1
            while max(hrf['range']) >> bits_req > 0:
                bits_req += 1
            hrf['bits'] = bits_req
            if hrf.has_key('default_value'):
                hrf['default_value'] = 0  # Default value if not defined
        return self.addrConf

    def _stas_parse_addrstr(self):
        '''
        Parse the address string
        '''
        assert isinstance(self.addrStr, str), "addrStr must be a string!"
        addrSpec = {}
        nBits = {}
        
        addrStrsplit = self.addrStr.split()
        addrStrsplit.reverse()
        
        for i in addrStrsplit:
            if i[0] in addrSpec:
                if len(addrSpec[i[0]]) < (int(i[1:]) + 1):
                    addrSpec[i[0]].extend( [None]*( int(i[1:])+1-len(addrSpec[i[0]]) ))
            else:
                addrSpec[i[0]] = [None] * (int(i[1:]) + 1)
                nBits[i[0]] = 0
                
            addrSpec[i[0]][int(i[1:])] = addrStrsplit.index(i)
        
        # Code for special ignore bits.
        # NOTE: Not necessary. Any letter can be used to ignore.
        nBitsTotal = 0
        for k, v in addrSpec.items():
            nBitsTotal += len(v)
            # Code for special ignore bits.
            # NOTE: Not necessary. Any letter can be used to ignore.
            #if k == 'I':
            #    addrSpec.pop(k)
            #    nBits.pop(k)
            #else:
            nBits[k] = len(v)
            
        return addrSpec, nBits, nBitsTotal

    def _stas_create_fields(self):
        """
        Creating fields for physical
        """
        self.nFields = len(self.nBits)
        self.field = [None for i in range(self.nFields)]
        for pos, pin in enumerate(self.addrPinConf):
            aspec = self.addrSpec[pin['id']]
            encoder = LayoutFieldEncoder(
                            aspec=np.array(aspec, 'uint'),
                            nWidth=self.nBits[pin['id']],
                            position=pos,
                            pin=pin['id'],
                            )
            self.field[pos] = encoder
        
        return self.field, self.nFields

