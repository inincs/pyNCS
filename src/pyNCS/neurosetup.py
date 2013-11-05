#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Sadique Sheik
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
from __future__ import with_statement
from xml.dom import minidom as md
from urllib2 import urlopen, URLError, HTTPError
from chip_v2 import NeuroChip
from mapping import Mapping, PMapping
from monitors import Monitors
import pyST
import warnings
from contextlib import contextmanager
from lxml import etree
from itertools import chain
from pyNCS.api import ComAPI
from pyNCS.api import ConfAPI


URL_SETUPDTD = 'http://ncs.ethz.ch/internal/files/setup.dtd/at_download/file'
URL_SETUPTYPEDTD = 'http://ncs.ethz.ch/internal/files/setuptype.dtd/at_download/file'
URL_TIMEOUT = 0.3
# temp


def dict_merge(x, y):
    '''
    Convenience function to merge two dictionaries
    returns the merged dictionary
    '''
    temp = x.copy()
    temp.update(y)
    return temp


def xml_parse_parameter(node):
    '''
    Convenience function to parse parameter leaves of XML
    '''
    d = dict()
    iterparams = node.iterfind('parameter')
    for i in iterparams:
        d[str(i.attrib['name'])] = i.text
    return d


def parse_and_validate(filename, dtd, validate = True):
    if 'http://' in filename:
        doc = etree.parse(urlopen(filename))
    else:
        doc = etree.parse(filename)
    
    if validate:
        try:
            parser = etree.DTD(urlopen(dtd, timeout=URL_TIMEOUT))
            parser.validate(doc)
            if len(parser.error_log.filter_from_errors()) > 0:
                for e in parser.error_log.filter_from_errors():
                    print e
                raise Exception('setupfile XML is not well formed, corrent the errors given above first')
            elif len(parser.error_log) > 0:
                for e in parser.error_log.filter_from_errors():
                    print e
        except URLError, e:
            warnings.warn('Could not access online setup.dtd (timeout {0}s), not validating'.format(URL_TIMEOUT))
        except HTTPError, e:
            warnings.warn('Could not access online setup.dtd (timeout {0}s), not validating'.format(URL_TIMEOUT))
        except Exception:
            warnings.warn('Unknown exception occured while validating setupfile (probably timed out)')

    return doc


def get_addrspec(node, typ):
    if node == None:
        return []
    else:
        if typ in node.attrib:
            return list(eval(node.attrib[typ]))
        else:
            return []


class NeuroSetup(object):
    def __init__(self,
            setuptype,
            setupfile,
            com_kwargs={},
            map_kwargs={},
            conf_kwargs={},
            prefix='./',
            offline=False,
            validate = True):
        '''
        Inputs:
        *setuptype*: Defines the channels and the slots in the setup. See NCS documentation.
        *setupfile*: Indicates which chips occupy the slots in the setup. See NCS documentation.
        *com_kwargs*: keyword arguments for communicator module. This gets merged with arguments in setup.xml
        *map_kwargs*: keyword arguments for communicator module. This gets merged with arguments in setup.xml
        *conf_kwargs*: keyword arguments for communicator module. This gets merged with arguments in setup.xml
        *prefix*: path prefix to setupfiles and chipfiles (default: current directory)
        *offline*: if True, the setup will not communicate, map or configure.
        *validate*: if True, the constructor will attempt to validate the setup and setuptype xml files.
        '''
        #Initialize and save API
        self.com_kwargs = com_kwargs
        self.conf_kwargs = conf_kwargs
        self.map_kwargs = map_kwargs
        #Save setup file locations
        self.setupfile = setupfile
        self.setuptype = setuptype
        self.prefix = prefix
        self.chips = {}
        self.chipslots = dict()
        self.slots = {}
        self.offline = offline
        self.validate = validate
        self.load_setuptype(self.setuptype, prefix = prefix, validate = validate)
        self.load(self.setupfile, prefix = prefix, offline = offline, validate = validate )
        self.aerDummyIn, self.aerDummyOut = self.aerDummy()
        self.update()
        self.apply()
        self.monitors = Monitors()

        if self.offline:
            warnings.warn('Running in Offline mode')

    def load_setuptype(self, filename, prefix='', validate = True):

        nsetup = parse_and_validate(filename, dtd=URL_SETUPTYPEDTD, validate = validate)

        for n in nsetup.iterfind('channelAddressing'):
#            if not n.getAttribute('name') == 'default':
#                print 'skipping channelAddressing %s'% n.getAttribute('name')
#                continue
            if str(n.attrib['type']) == 'monitor':
                self.monBits = eval(n.attrib['bits'])
            elif str(n.attrib['type']) == 'sequencer':
                self.seqBits = eval(n.attrib['bits'])
            else:
                raise TypeError('channelAddressing type should be either monitor or sequencer')

        self.aerSlot = dict()
        for nslot in nsetup.iterfind('slot'):
            #Considfer defining a function (for esthetic reasons)
            id = int(nslot.attrib['id'])
            self.aerSlot[id] = dict()
            ######### Mon
            self.aerSlot[id]['monIn'] = get_addrspec(nslot.
                find('aerMon'), 'in')
            self.aerSlot[id]['monOut'] = get_addrspec(nslot.
                find('aerMon'), 'out')

            ######### Seq
            self.aerSlot[id]['seqIn'] = get_addrspec(nslot.
                find('aerSeq'), 'in')
            self.aerSlot[id]['seqOut'] = get_addrspec(nslot.
                find('aerSeq'), 'out')

    def load(self, filename, prefix='', offline=False, validate = True):
        '''
        Loads the setup
        Inputs:
        *filename*: setup file name
        *prefix*: path to be prepended to chipfile names
        *offline*: if True, the chips will not be configured ("pretend" mode).
        '''
        nsetup = parse_and_validate(filename, dtd=URL_SETUPDTD, validate = validate)

        #parse defaultchip (should be unique)
        self.defaultchipfile = str(
            nsetup.find('defaultchip').attrib['chipfile'])
        #
        # Load communicator, currently, only only communicator per setup is
        # supported
        for ncom in nsetup.iterfind('communicator'):
            com_kwargs = xml_parse_parameter(ncom)
            self.com_kwargs = dict_merge(self.com_kwargs, com_kwargs)
            try:
                self.com_api = self._import_module(str(ncom.attrib['module']))
            except ImportError, e:
                warnings.warn(e.message)
                self.com_api = ComAPI
                self.com_kwargs = {}
            self.communicator = self.com_api.Communicator(**self.com_kwargs)
        
        # Load virtual chips
        for nchip in nsetup.iterfind('virtualchip'):
            chipid = str(nchip.attrib['id'])
            chipfile = prefix + str(nchip.attrib['chipfile'])
            slot = int(eval(nchip.attrib['slot']))
            chip = NeuroChip(chipfile, id=chipid, offline=True)
            #Be sure that the chip is virtual
            chip.virtual = True
            self.chips[chipid] = chip
            self.chipslots[chipid] = slot
            self.slots[slot] = chipid

        #Load configurators and chips
        for nchip in nsetup.iterfind('chip'):
            chipid = str(nchip.attrib['id'])
            slot = int(eval(nchip.attrib['slot']))
            chipfile = prefix + str(nchip.attrib['chipfile'])
            nconf = nchip.find('configurator')
            module = str(nconf.attrib['module'])
            conf_kwargs = xml_parse_parameter(nconf)
            try:
                conf_api = self._import_module(module)
            except ImportError, e:
                warnings.warn(e.message)
                conf_api = ConfAPI
                conf_kwargs = {}
            chip = NeuroChip(chipfile, id=chipid, offline=offline,
                             conf_api=conf_api,
                             conf_kwargs=conf_kwargs)
            self.chips[chipid] = chip
            self.chipslots[chipid] = slot
            self.slots[slot] = chipid
            chip.configurator.register_neurosetup(self)

        #Load Mapper
        # load a virtual mapping table by defauld
        self.mapping = Mapping('Virtual Mapping') 

        for nmapper in nsetup.iterfind('mapper'):
            map_kwargs = xml_parse_parameter(nmapper)
            self.map_kwargs = dict_merge(self.map_kwargs, map_kwargs)
            try:
                self.map_api = self._import_module(str(nmapper.attrib['module']))
            except ImportError, e:
                warnings.warn(e.message)
                self.map_api = ConfAPI
                self.map_kwargs = {}
            self.mapper = self.map_api.Mappings(**self.map_kwargs)
            if self.map_kwargs.has_key('version'):
                if float(self.map_kwargs['version']) == 3.0:
                    self.mapping = PMapping('mapping')
                else:
                    self.mapping = Mapping('mapping')

        #If there is no mapper tag in setup, build an empty mapper
        if not hasattr(self, 'mapper'):
            from ConfAPI import Mappings
            self.mapper = Mappings()

    def _import_module(self, module):
        try:
            #without fromlist, only the package is returned
            mod = __import__(module, fromlist=[module])
        except ImportError, error:
            #print(error)
            raise ImportError('{0} API failed to be imported'.format(module))
        return mod

    def aerDummy(self):
        ''' returns a placeholder pyST.addrSpec '''
        if self.defaultchipfile.endswith('.csv'):
            aerIn, aerOut = pyST.STas.load_stas_from_csv(self.prefix+\
                                                     self.defaultchipfile)
        elif self.defaultchipfile.endswith('.nhml'):
            aerIn, aerOut = pyST.STas.load_stas_from_nhml(self.prefix+\
                                                     self.defaultchipfile)
        return aerIn, aerOut

    def update(self):
        '''
        updates the default monitor/sequencer (pyST) with the chips contained
        in this setup -- always call this function when finished with adding chips
        '''
        self.chipsIn = [None for i in range(2 ** len(self.seqBits))]
        self.chipsOut = [None for i in range(2 ** len(self.monBits))]
        seqList = [self.aerDummyIn for i in range(2 ** len(self.seqBits))]
        monList = [self.aerDummyOut for i in range(2 ** len(self.monBits))]

        for id, nslot in self.chipslots.iteritems():
            for ii in self.aerSlot[nslot]['monIn']:
                try:
                    if self.chips[id].aerIn is not None:
                        if monList[ii] is self.aerDummyOut or \
                           not self.chips[id].virtual:
                            monList[ii] = self.chips[id].aerIn
                            self.chipsIn[ii] = self.chips[id]
                except KeyError, e:
                    pass
                #If mon.aerIn doesn't exist, then skip it
                except AttributeError, e:
                    pass

            for ii in self.aerSlot[nslot]['monOut']:
                try:
                    if self.chips[id].aerOut is not None:
                        if monList[ii] is self.aerDummyOut or \
                           not self.chips[id].virtual:
                            monList[ii] = self.chips[id].aerOut
                            self.chipsOut[ii] = self.chips[id]
                except KeyError, e:
                    raise e
                except AttributeError, e:
                    pass

            for ii in self.aerSlot[nslot]['seqIn']:
                try:
                    if self.chips[id].aerIn is not None:
                        if seqList[ii] is self.aerDummyIn or \
                           not self.chips[id].virtual:
                            seqList[ii] = self.chips[id].aerIn
                            self.chipsIn[ii] = self.chips[id]
                except KeyError, e:
                    raise e
                except AttributeError, e:
                    pass
                #If seq.aerIn doesn't exist, then skip it

            for ii in self.aerSlot[nslot]['seqOut']:
                try:
                    if self.chips[id].aerOut is not None:
                        if seqList[ii] is self.aerDummyIn or \
                           not self.chips[id].virtual:
                            seqList[ii] = self.chips[id].aerOut
                            self.chipsOut[ii] = self.chips[id]
                except KeyError, e:
                    raise e
                except AttributeError, e:
                    pass

        self.seq = pyST.channelAddressing(channelBits=self.seqBits, stasList=seqList)
        self.mon = pyST.channelAddressing(channelBits=self.monBits, stasList=monList)

    def apply(self):
        ''' sets default monitor/sequencer to this setup '''
        #?? self.update()
        pyST.setDefaultMonChannelAddress(self.mon)
        pyST.setDefaultSeqChannelAddress(self.seq)

    def prepare(self):
        self.mapping.prepare()
        if not self.offline:
            if len(self.mapping.mapping) > 0:
                self.mapper.set_mappings(self.mapping.mapping)

    def __copy__(self):
        return self.__class__(self.setuptype, self.setupfile,
                              prefix=self.prefix, offline=self.offline)

    def __deepcopy__(self, memo):
        return self.__copy__()

    def __getstate__(self):
        return {'setuptype': self.setuptype,
                "setupfile": self.setupfile,
                "prefix": self.prefix,
                "offline":self.offline,
                "validate": self.validate,
               }

    def reload(self):
        '''
        Call this function if you recovered the object from a pickle and doesnt
        have any of the information or the information is corrupted.
        '''
        self.com_kwargs = {}
        self.conf_kwargs = {}
        self.map_kwargs = {}
        self.chips = {}
        self.chipslots = {}
        self.slots = {}
        self.load_setuptype(self.setuptype, prefix=self.prefix, validate = self.validate)
        self.load(self.setupfile, prefix=self.prefix, offline = self.offline, validate = self.validate)
        self.aerDummyIn, self.aerDummyOut = self.aerDummy()
        self.update()

    def __getinitargs__(self):
        return self.setuptype, self.setupfile

    def _pre_process(self, stim):
        evs_in = self.seq.exportAER(stim, isi=True)
        return evs_in.get_adtmev()

    def _post_process(self, evs, filter_channels=None):
        evs_out = pyST.events(evs, 'p')
        mon_ch_addr = self.mon
        #extract per channel events -> ch_events
        ch_evs = mon_ch_addr.extract(evs_out)
        #Filter
        #-Channel Events filtering is easier to use than RawOutput.decode_dict
        if filter_channels != None:
            ch_evs.filter_channel(filter_channels)
        stout = mon_ch_addr.rawoutput_from_chevents(
            ch_evs,
            normalize=False,
            filter_duplicates=False
            )
        self.monitors.populate_monitors(chstlist=stout)
        return stout

    def run(self, *args, **kwargs):
        '''
        Prepares and stimulates.
        *args* and *kwargs* are keyword arguments passed to Communicator.run()
        '''
        self.prepare()
        if not self.offline:
            return self.stimulate(*args, **kwargs)

    def stimulate_raw(self, raw_stim, **kwargs):
        '''
        Calls communicator.run without pre- or post-processing
        *raw_stim*: a numpy array in (addr, time) format. Type should be uint32 (shape = (-1,2)).
        Useful for debugging purposes. Format of raw_stim corresponds to the one defined in ComAPI:

        '''
        return self.communicator.run(raw_stim, **kwargs)

    def stimulate(self, stim={}, **kwargs):
        '''
        Run without preparing.
        Pre-processes, runs communicator, and post-processes.
        *kwargs* are keyword arguments passed to self.communicator.run()
        Returns a Stas.RawOutput object and populates monitors. The latter is the preferred way of reading data out.
        '''
        stim_evs = self._pre_process(stim)
        #run rec for run and record (consider using run for recording TODO)
        evs = self.communicator.run_rec(stim_evs, **kwargs)
        return self._post_process(evs, self.monitors.channels)


