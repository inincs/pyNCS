#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Emre Neftci
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
# -*- coding: utf-8 -*-
######################################################################
# Emre Neftci Author <emre.neftci@ini.phys.ethz.ch>
#######################################################################
from __future__ import absolute_import
import numpy as np
import os
import copy
import time
import warnings
from .STsl import *
import itertools
from contextlib import contextmanager
from . import pyST_globals
from lxml import etree


PERFORMANCE_DEBUG=False

#TODO: RawOutput should return spikelists with complete id_list even when empy
class RawOutput(object):
    '''
    RawOutPut is a class which contains raw AER data per channel (Physical addresses).

    It decodes the AER data "just in time"

    Inputs:
    *raw_data*: events with physical addresses
    *id_data*: all addresses in given channel
    *decoder_dict*: Dictionary of functions for decoding addresses. Key designates the channel
    *t_start*: start time of all spikelists in ms
    *t_stop*: end time of all spikelists in ms
    *filter_duplicates*: If True erases in all channels double events within the a 0.01ms time-frame (buggy hw)

    Usage:
    #In conjunction with pyAex.netClient

    >>> client = client.netClient(host='localhost')
    >>> raw_out = client.stimulate()
    >>> client.stop()
    >>> for sl in raw_output: sl.raster_plot()
    '''
    def __init__(self, raw_data, decoder_dict, t_start=0, t_stop=1000, filter_duplicates=False):

        #Inputs
        self.raw_data = raw_data.copy()
        self.decoder_dict = decoder_dict
        self.filter_duplicates = filter_duplicates

        #Containers
        self.decoded_data = {}
        self.channels = self.raw_data.keys()
        self.t_start = t_start
        self.t_stop = t_stop
        self.print_num_events()

    def print_num_events(self):
        lev = []
        for k, v in self.raw_data.iteritems():
            lev.append('Ch{0}: {1} evs '.format(k, len(v)))
        print(''.join(lev))

    def __getitem__(self, key):
        key = int(key)
        if not key in self.decoded_data:
            self.decode_data(key)

        return self.decoded_data[key]

    def __getstate__(self):
        self.decode_all_channels()
        dc = self.__dict__.copy()
        del dc['decoder_dict']
        return dc

    def decode_all_channels(self):
        for i in self.channels:
            self.decode_data(i)

    def decode_data(self, key):
        with self.check_has_key_somewhere(key):
            if PERFORMANCE_DEBUG: t0=time.time()
            ad_data = self.decoder_dict[key](self.raw_data[key].get_ad())
            tm_data = self.raw_data[key].get_tm().astype('float') / 1000
            if PERFORMANCE_DEBUG:
                print('Decoding events took {0} seconds'.format(time.time()-t0))

            if PERFORMANCE_DEBUG: t0=time.time()
            evs = events(atype='l')
            evs.add_adtm(ad_data, tm_data)
            if PERFORMANCE_DEBUG:
                print('Building events took {0} seconds'.format(time.time()-t0))

            if PERFORMANCE_DEBUG: t0=time.time()
            st_data = SpikeList( evs.get_adtmev(), np.unique(evs.get_ad()))
            if PERFORMANCE_DEBUG:
                print('Building SpikeList took {0} seconds'.format(time.time()-t0))

            if self.filter_duplicates:
                st_data.filter_duplicates()
            self.raw_data.pop(key)
            self.decoded_data[key] = st_data

    def iteritems(self):
        for key, value in self.iterchannels(self.channels):
            yield key, value

    def __iter__(self):
        for key, value in self.iterchannels(self.channels):
            yield value

    def iterchannels(self, channels):
        for key in channels:
            #key is a channel
            #self[key] is a spikelist
            yield key, self[key]

    @contextmanager
    def check_has_key_somewhere(self, key):
        try:
            yield
        except KeyError:            
            if not key in self.decoder_dict:
                raise KeyError("There is no function to decode %d" % key)
            if key in self.decoded_data:
                print("Channel {0} has already been decoded".format(key))
            elif not key in self.raw_data:
                #print("Channel {0} is not present, assuming no events".format(key
                self.decoded_data[key] = SpikeList([], [])
            else:
                raise


class events(object):
    #Number of dimensions (address, timestamp)
    NF = 2

    def __init__(self, ev=None, atype='p', isISI=False):
        self.isISI = isISI
        assert isinstance(atype, str)
        __atype = atype[0].lower()

        if __atype == 'p':
            __dtype = np.dtype([('tm', 'uint32'),
                              ('ad', 'uint32')])

        if __atype == 'l':
            __dtype = np.dtype([('tm', 'float'),
                              ('ad', 'float')])

        self.__atype = __atype
        self.__dtype = __dtype

        if isinstance(ev, events):
            self.__data = ev.__data.copy()
        elif ev is not None:
            ev = np.array(ev)
            if ev.shape[1] == self.NF:
                self.set_data(ev[:, 0], ev[:, 1])
            elif ev.shape[0] == self.NF:
                self.set_data(ev[0, :], ev[1, :])
            else:
                raise TypeError
        else:
            self.__data = np.zeros([0], self.dtype)

    @property
    def data(self):
        return self.__data

    @property
    def atype(self):
        return self.__atype

    @property
    def dtype(self):
        return self.__dtype

    def __len__(self):
        return self.get_nev()

    def __add__(self, other):
        self.add_adtm(other.ad, other.tm)

    def __repr__(self):
        return self.get_adtm().__repr__()

    def get_nev(self):
        return len(self.get_tm())

    @property
    def nev(self):
        return self.get_nev()

    def get_tdur(self):
        if self.nev == 0:
            return 0
        if self.isISI:
            return self.tm.sum()
        else:
            if self.get_nev() > 0:
                return self.tm[-1] - self.tm[0]

    @property
    def tdur(self):
        return self.get_tdur()

    def get_t_stop(self):
        if self.isISI:
            return self.tm.sum()
        else:
            if self.get_nev() > 0:
                return self.tm[-1]

    @property
    def t_stop(self):
        return self.get_t_stop()

    def get_ad(self):
        return self.__data['ad']

    def set_ad(self, ad):
        self.__data['ad'] = ad
        return self.get_ad()

    @property
    def ad(self):
        return self.get_ad()

    #@ad.setter
    #def ad(self, value):
    #    return self.set_ad(value)

    def get_tm(self):
        return self.__data['tm']

    def set_tm(self, tm):
        self.__data['tm'] = tm
        return self.get_tm()

    @property
    def tm(self):
        return self.get_tm()

    #@tm.setter
    #def tm(self, value):
    #    return self.set_tm(value)

    def set_data(self, ad, tm):
        assert len(ad) == len(tm), "addresses and timestamps lengths are incompatible %d %d" % (len(ad), len(tm))
        self.__data = np.zeros(len(ad), self.dtype)
        self.__data['ad'] = ad.astype(self.dtype['ad'])
        self.__data['tm'] = tm.astype(self.dtype['tm'])

    def add_adtmev(self, ev):
        if not isinstance(ev, np.ndarray):
            ev = np.array(ev)

        if len(ev.shape) != self.NF:
            ev = ev.reshape(-1, 2).astype(self.dtype)

        ad = np.concatenate([self.ad, ev[:, 0]])
        tm = np.concatenate([self.tm, ev[:, 1]])
        self.set_data(ad, tm)

    def add_adtm(self, ad, tm):
        if not isinstance(ad, np.ndarray):
            ad = np.array(ad)

        if len(ad.shape) != 1:
            ad = ad.reshape(-1, 1)

        if not isinstance(tm, np.ndarray):
            tm = np.array(tm)

        if len(tm.shape) != 1:
            tm = tm.reshape(-1, 1)

        assert tm.shape == ad.shape

        ad = np.concatenate([self.ad, ad])
        tm = np.concatenate([self.tm, tm])
        self.set_data(ad, tm)

    def get_tmad(self):
        return np.array([self.tm, self.ad])

    def get_adtm(self):
        return np.array([self.ad, self.tm])

    def get_tmadev(self):
        return self.get_tmad().transpose()

    def get_adtmev(self):
        return self.get_adtm().transpose()

    def normalize_tm(self, t0=0.):
        if not self.isISI:
            t_start = self.get_tm().min()
            self.set_tm(self.get_tm() - t_start + t0)
        else:
            t = self.get_tm()
            t[0] = t0
            self.set_tm(t)

    def get_adisi(self):

        if self.isISI:
            return self.get_adtm()
        else:
            if self.nev > 0:
                tm = np.concatenate([np.array([self.tm[0]]), np.diff(self.tm)])
                return np.array([self.ad, tm])
            else:
                return np.zeros([2,0])

    def get_isiad(self):

        if self.isISI:
            return self.get_tmad()
        else:
            if self.nev > 0:
                tm = np.concatenate([np.array([self.tm[0]]), np.diff(self.tm)])
                return np.array([tm, self.ad])
            else:
                return np.zeros([2,0])

    def set_isi(self):
        if self.isISI:
            pass
        else:
            evs = self.get_adisi()
            self.set_data(evs[0],evs[1])
            self.isISI = True
            

    def set_abs_tm(self):
        '''
        Transform ISI timestamps into absolute time 
        '''
        if self.isISI:
            self.ad, self.set_tm(np.cumsum(self.get_tm()))
            self.isISI = False
        else:
            pass

    def filter_by_mapping(self, mapping):
        """
        Map the events, given a mapping dictionary like:
        map[src]=[target1,target2,...,targetn],
        """
        #The following is optimized for performance
        wasISI = False
        if self.isISI:
            wasISI = True
            self.set_abs_tm()
        evs = self.get_adtmev()
        ff = lambda x: x[0] in mapping
        filt = filter(ff, evs) #keep only addresses that are mapped
        if len(filt) > 0:
            evs_filt = events(filt, self.atype).get_adtmev()
            #Get mapped addresses
            #list(chain(* concatenates lists of lists for low cost (i.e. O(n))
            m_ad = np.array(list(itertools.chain(*map(mapping.get, evs_filt[:, 0]))), self.dtype['ad'])
            #Get the respective timestamps
            m_tm = np.array(list(itertools.chain(*map(lambda x: len(mapping.get(x[0])) * [x[1]], evs_filt))), self.dtype['tm'])

            self.set_data(m_ad,m_tm)
        else:
            self.empty()
        if wasISI:
            self.set_isi()

    def empty(self):
        self.__data = np.zeros([0], self.dtype)

    def iter_by_timeslice(self, tm):
        import bisect
        #ISI -> Cumulative
        if self.isISI:
            sum_tm = np.cumsum(self.get_tm())
        else:
            sum_tm = self.get_tm()
        #No in-place change
        id_start = 0
        t = 0

        #better to recycle an object rather than creating new (slow)
        evs = events(isISI=self.isISI)

        while id_start < len(sum_tm):
            t += tm
            id_stop = bisect.bisect_right(sum_tm, t, lo=id_start)
            evs.__data = self.__data[id_start:id_stop]
            id_start = id_stop
            rest = tm - evs.get_tdur()
            print(tm, evs.get_tdur())

            if not evs.get_nev() > 0:
                rest = 0

            t -= rest
            #ISI or not is determined
            yield evs, tm - rest

    def __iter__(self):
        for ev in self.get_tmadev():
            yield ev[0],ev[1]

    def sort(self):
        '''
        Sort events by addresses and timestamps (in this order).
        '''
        self.set_abs_tm()
        self.__data.sort(order=['ad', 'tm'])

    def sort_tm(self):
        '''
        Sort events by timestamps and addresses (in this order).
        '''
        self.set_abs_tm()
        self.__data.sort(order=['tm', 'ad'])

    def demultiplex(self):
        '''
        Generates a dictionary with addesses as keys and a list of timestamps as values.
        Used internally for generating SpikeLists
        '''
        evs = events(ev=self)
        evs.sort()
        ads = np.unique(evs.ad)

        d = dict()
        k_start = 0
        k_stop = np.searchsorted(self.ad, ads, side='right')

        for i, a in enumerate(ads):
            d[a] = self.tm[k_start:k_stop[i]]
            k_start = k_stop[i]

        return d


class channelEvents(dict):
    '''
    inputs:
    *channel_events*: dictionary or channelEvents.
    *atype*: Address type 'physical' or 'logical'
    '''
    def __init__(self, channel_events=None, atype='physical'):
        assert isinstance(atype, str)
        self.__atype = atype.lower()[0]

        if isinstance(channel_events, dict):
            for k, v in channel_events.iteritems():
                self.add_ch(k, events(v, self.atype))
        elif channel_events is not None:
            raise TypeError("channel_events must be a dictionary, None or a channelEvents. Alternatively, use an events object and extract from channelAddressing")

    @property
    def atype(self):
        return self.__atype

    def copy(self):
        return channelEvents(self, atype=self.atype)

    def __add__(self, other):
        for i in other.keys():
            if i in self:
                self[i] + other[i]
            else:
                self[i] = other[i]

    def __getattr__(self, attrName):
        """
        Check if the attribute exists in self, otherwise look for attributes of the encapsulated events and build a function by adding a channel argument to it.
        """
        if attrName in self.__dict__:
            return self.__dict__.__getitem__(attrName)
        else:
            def dnfun(ch, *args, **kwargs):
                return getattr(self[ch], attrName)(*args, **kwargs)
            return dnfun

    def add_ch(self, channel, ev):
        '''
        Add an events object to a channel
        Remark: Creates a reference to the provided events! In place changes also affect the events!
        '''
        assert channel not in self
        if not isinstance(ev, events):
            self[channel] = events(ev, atype=self.atype)
        else:
            self[channel] = ev

    def __len__(self):
        return self.get_nev()

    def add_adtmch(self, channel, ad, tm):
        if channel not in self:
            self[channel] = events(atype=self.atype)

        self.add_adtm(channel, ad, tm)

    def flatten(self):
        ev = events(atype=self.atype)
        for ch in self:
            ev.add_adtmev(self.get_adtmev(ch))
        return ev

    def filter_channel(self, channel_list=None):
        """
        Removes (in-place) all channels with are not in channel_list.
        If channel_list is omitted, the ch_events is left unchanged
        """
        if channel_list == None:
            return None
        all_ch = self.keys()
        for ch in all_ch:
            if ch not in channel_list:
                self.pop(ch)

    def get_all_tm(self):
        return self.flatten().get_tm()

    def get_last_tm(self):
        t = 0
        for evs in self.itervalues():
            t=max(t,evs.get_tm()[-1])
        return t

    def get_first_tm(self):
        t = np.inf
        for evs in self.itervalues():
            t=min(t,evs.get_tm()[0])
        return t

    def get_all_ad(self):
        return self.flatten().get_ad()

    def get_all_adtmev(self):
        return self.flatten().get_adtmev()

    def get_all_tmadev(self):
        return self.flatten().get_tmadev()

    def get_nev(self):
        n = 0
        for v in self.itervalues():
            n+=len(v)
        return n

    def iter_by_timeslice(self, tm):
        return self.flatten().iter_by_timeslice(tm)

    def filter_all_by_mapping(self, mapping):
        """
        Modifies, in-place, the encapsulated events acccording to the one-to-many mapping (key is int/float, value is iterable)
        """
        for ch in self:
            self[ch].filter_by_mapping(mapping)

    def filter_all_by_channel_mapping(self, mapping):
        """
        Modifies, in-place, the encapsulated events acccording to the one-to-many mapping (key is int/float, value is iterable)
        In this function, a different mapping is used for each channel (useful for address that do not contain channel information, for example)
        """
        for ch in self:
            self[ch].filter_by_mapping(mapping[ch])


def setDefaultMonChannelAddress(cs):
    """ Sets the default Monitoring Channel Addressing scheme (i.e. AER Output) in pyST globals, and used by pyAex. The argument is then returned when getDefaultMonChannelAddress is called

    cs: The Channel Addressing Scheme, instance of the channelAddressing class

    See also:
    `getDefaultMonChannelAddress <#pyST.STas.getDefaultMonChannelAddress>`_
    """

    pyST_globals.DefaultMonChannelAddress = cs
    return None


def getDefaultMonChannelAddress():
    """
    Returns the default Monitoring Channel Addressing Scheme (i.e. AER Output). No arguments

    See also:
    `setDefaultMonChannelAddress <#pyST.STas.setDefaultMonChannelAddress>`_
    """

    try:
        pyST_globals.DefaultMonChannelAddress
    except NameError:
        print('Error, run setDefaultMonChannelAddress first')
    return pyST_globals.DefaultMonChannelAddress


def setDefaultSeqChannelAddress(cs):
    """
    Sets the default Sequencer Channel Addressing scheme (i.e. AER Input) in pyST globals, and used by pyAex. The argument is then returned when getDefaultMonChannelAddress is called

    cs: The Channel Addressing Scheme, instance of the channelAddressing class

    See also:
    `pyST.getDefaultSeqChannelAddress <#pyST.STas.getDefaultSeqChannelAddress>`_
    """

    pyST_globals.DefaultSeqChannelAddress = cs
    return None


def getDefaultSeqChannelAddress():
    """
    Returns the default Monitoring Channel Addressing Scheme (i.e. AER Input). No arguments

    See also:
    `setDefaultSeqChannelAddress <#pyST.STas.setDefaultSeqChannelAddress>`_
    """

    try:
        pyST_globals.DefaultSeqChannelAddress
    except NameError:
        print('Error, run setDefaultSeqChannelAddress first')
    return pyST_globals.DefaultSeqChannelAddress


# Field class used for addrSpecIgnoreSynapseNeuron and STIgnoreNeuronSynapse
class channelAddressing:
    """
    This class contains the information defining the Channel Addressing Scheme for a multichip setup, i.e. definition of which chips are found on which channel

    channelAddressing(stasList,nChannelBits)

    *nChannelBits*: Positions of channel bits in the hardware addresses. Default nChannelBits=[15,16]

    *stasList*: list of address specification objects, of the same length as len(channelBits)**2

    See also:
    `addrSpec <#pyST.STas.addrSpec>`_

    """

    def __init__(self, stasList, channelBits=[15, 16]):
        """
        Constructor of the channelAddressing object

        See also
            SpikeTrain
        """
        self.nBitsTotal = np.array([0] * len(stasList), 'uint16')
        self.nChannelBits = len(channelBits)
        self.nChannels = 2 ** len(channelBits)
        self.channels = range(self.nChannels)

        if len(stasList) <= self.nChannels:
            for i in xrange(len(stasList)):
                self.nBitsTotal[i] = stasList[i].nBitsTotal
            for i in xrange(len(stasList), self.nChannels):
                stasList[i].append(None)
        else:
            raise RuntimeError("len(stasList)>2**nBitChannel !")

        self.nBitsTotal = min(channelBits) - 1
        self.stasList = stasList

    def __getitem__(self, channel):
        return self.stasList[channel]

    def __len__(self):
        """
        Returns the length of the channel addressing scheme
        """

        return len(self.stasList)

    def __repr__(self):
        repr_str = "\n"
        for i in range(len(self)):
            repr_str += "Channel " + str(i) + ": " + self[i].__repr__() + "\n"
        return repr_str

    def __iter__(self):
        for addr_spec in self.stasList:
            yield addr_spec

    def reprAddressSpecification(self):
        '''
        String representation of the channel address specification
        '''
        base_label  = '|--Channel--|----Chip-----|\n'
        base_string = '|--{0}bits----|----{1}bits---|\n\n'.format(self.nChannelBits,self.nBitsTotal)
        chann_strings = []
        for n, addr_spec in enumerate(self):
            chann_strings.append('Channel {0}:\n'.format(n))
            chann_strings.append(addr_spec.repr_addr_spec(self.nBitsTotal))
            chann_strings.append('\n')
        repr_addr_spec = base_label+base_string+''.join(chann_strings)
        return repr_addr_spec


    def getValue(self, channel):
        """
        returns the channel mask

        channel: int between 0 and len(channelBits)**2
        """

        return np.uint32(channel << self.nBitsTotal)

    def extract(self, ev):
        """
        Extracts the channel information from an array of hardware events. Returns a channelEvents object. Use addrPhysicalExtract for decoding the adresses instead.

        *eventsArray*: numpy array of hardware events (uint32)

        See also:
            addrPhysicalExtract
        """
        ch_events = channelEvents(atype='Physical')
        channel = ev.get_ad() >> self.nBitsTotal     # Get channel information.

        for channelIdx in xrange(self.nChannels):
            t = pylab.find(channelIdx == channel)
                 # Much faster than boolean list or filter
            if len(t) > 0:
                ad = (ev.get_ad()[t]) & (2 ** self[channelIdx].nBitsTotal - 1)
                ch_events.add_adtmch(channelIdx, ad, ev.get_tm()[t])
        return ch_events

    def isChannelAddrList(self, addr):
        """
        Mostly internal function


        Check whether addr is a valid channelAddr list/dict

        *addr*: Either list of dimension channelAddressing.nChannels, containing arrays of addresses. The list index is the channel number, or dict whose keys are the channel numbers.

        See also:

            Stas.isValidAddress
        """
        if isinstance(addr, list):
            assert len(addr) == self.nChannels
        elif isinstance(addr, dict):
            possible_channels = range(self.nChannels)
            for i in addr.iterkeys():
                assert i in possible_channels, "The keys of the addr dict mst be valid channel numbers, i.e. between %d and %d" % (0, self.nChannels)

            addr_list = [None for i in xrange(self.nChannels)]
            for k, a in addr.iteritems():
                addr_list[k] = a
            addr = addr_list
        else:
            raise RuntimeError("addr must be a list or a dict")
        return addr

    def addrLogicalConstruct(self, addr):
        """
        Constructs Logical addresses, *i.e.* addresses in float format useful for plotting because they keep the neuron - synapse order.

        Returns a float numpy array.

        **NOTE:** Logical addresses do not contain channel information

        *addr*: address (in a form isChannelAddrList() can understand), such as ``{0:[range(15),2]}``

        See also:

            addrLogicalConstruct, addrLogicalExtract
        """
        addr = self.isChannelAddrList(addr)
        mainAddr = [None for i in xrange(self.nChannels)]

        #Construct channel by channel if not empty
        for channelIdx in range(len(addr)):
            if addr[channelIdx] is not None:
                mainAddr[channelIdx] = np.array(self[channelIdx]
                    .addrLogicalConstruct(addr[channelIdx]), 'float')

        return mainAddr

    def addrLogicalExtract(self, addr):
        """
        Extracts Logical Address with channel information

        *addr*: address (in a form isChannelAddrList() can understand), such as ``{0:[1.,2.,3.]}``
        """
        addr = self.isChannelAddrList(addr)
        mainAddr = [None for i in xrange(self.nChannels)]

        #Construct channel by channel if not empty
        for channelIdx in xrange(self.nChannels):
            if addr[channelIdx] is not None:
                addr[channelIdx] = np.array(addr[channelIdx], 'float')
                mainAddr[channelIdx] = self[
                    channelIdx].addrLogicalExtract(addr[channelIdx])

        return mainAddr

    def addrPhysicalConstruct(self, addr):
        """
        Constructs Physical addresses to human readable addresses

        *addr*: dictionary of human readable numbers ( the physical addresses ), with channel numbers as keys

        """
        addr = self.isChannelAddrList(addr)
        mainAddr = np.zeros([0], dtype='uint32')
        for channelIdx in xrange(self.nChannels):
            if addr[channelIdx] is not None:
                mainAddr = np.concatenate((
                        mainAddr,
                        self[channelIdx].addrPhysicalConstruct(
                            addr[channelIdx])
                        + self.getValue(channelIdx
                        )))
        #Add Channel Part
        return mainAddr

    def addrPhysicalExtract(self, addr):
        """
        Extracts Physical addresses to human readable addresses

        *addr*: numpy array of uint32 numbers ( the physical addresses )

        Output is channelAddrList, channel where channelEventsList is a list of arrays which contain the addresses and channel is an array of channels.
        """
        if not isinstance(addr, np.ndarray):
            addr = np.array(addr, np.uint32)
        if not addr.dtype == np.uint32:
            addr = addr.astype(np.uint32)

        channels_in_addr = addr >> self.nBitsTotal
        channelEventsList = [None for i in xrange(self.nChannels)]
        for channelIdx in np.unique(channels_in_addr):
            t = (channels_in_addr == channelIdx)
            channelEventsList[channelIdx] =\
                    self[channelIdx].addrPhysicalExtract(
                        addr[t] - (channels_in_addr[t] << self.nBitsTotal))
        return channelEventsList

    def addrPhysicalLogical(self, addr):
        """
        Extracts Physical addresses to logical addresses (directly)

        *addr*: numpy array of uint32 numbers ( the physical addresses )

        Output is channelAddrList, channel where channelEventsList is a list of arrays which contain the addresses and channel is an array of channels.
        """
        if not isinstance(addr, np.ndarray):
            addr = np.array(addr, np.uint32)
        if not addr.dtype == np.uint32:
            addr = addr.astype(np.uint32)

        channels_in_addr = addr >> self.nBitsTotal
        channelEventsList = [None for i in xrange(self.nChannels)]
        for channelIdx in np.unique(channels_in_addr):
            t = (channels_in_addr == channelIdx)
            channelEventsList[channelIdx] =\
            self[channelIdx].addrPhysicalLogical(addr[t] -
                (channels_in_addr[t] << self.nBitsTotal))
        return channelEventsList

    def importAER(self, input=None, sep='\t', dt=1e-6, format='a', isi=False, *args, **kwargs):
        """
        Function for extracting, translating events from a numpy array. Output is a channelEvents object.
        
        Inputs:
        *input*: if a string, will be treated as a filename and passed to np.loadtxt, if a numpy array, will be considered as events of dimension 2 x number of events. By default, the addresses are on [:,0]
        *format*: either 't' or 'a' respectively meaning timestamps and addresses on the first column. ('a' by default)
        *kwargs*: keyword arguments passed to np.loadtxt
        """

        if isinstance(input, str):
            try:
                if 'converters' in kwargs:
                    ae = np.loadtxt(input, **kwargs)
                else:
                    ae = np.loadtxt(
                        input, delimiter=sep, dtype='uint32', **kwargs)
            except IOError:
                ae = np.zeros([0, 2], 'uint32')

            if format == 'a':
                input = events(ae, 'p')
            elif format == 't':
                input = events(np.fliplr(ae), atype='Physical')
            else:
                raise RuntimeError("format must be a (addresses first) or t (timestamps first)")
        elif isinstance(input, channelEvents):
            pass
        else:
            input = events(input)

        if isi == True:
            input.set_tm(np.cumsum(input.get_tm()))

        ch_events = self.extract(input)
        ch_events_log = channelEvents(atype='Logical')

        for ch in ch_events.iterkeys():
            if len(ch_events.get_tm(ch)) > 0:
                ad = self[ch].addrPhysicalLogical(ch_events.get_ad(ch))
                ch_events_log.add_adtmch(ch, ad, ch_events.get_tm(ch))

        return ch_events_log

    def normalizeAER(self, ch_events):
        """
        Called before extract to throw away pre-stimulus data
        """
        if not isinstance(ch_events, channelEvents):
            raise TypeError("ch_events must be a channelEvents object")

        # Find the starting time. Empty spike trains are accepted with None and
        # shape[0]=0 arrays
        #Normalize
        if ch_events.get_nev() == 0:
            raise RuntimeError("ch_events does not contain any events!")

        tStartGlobal = ch_events.get_all_tm().min()
        for ch in ch_events:
            ch_events.set_tm(ch, ch_events.get_tm(ch) - tStartGlobal)
        return ch_events

    def generateST(self, ch_events, normalize=True):
        """
        Extracts events from eventsChannel, using the address specification specified for the given channels, returns a list of SpikeList
        """
        assert ch_events.atype == 'l', 'channelEvents must be of type logical'

        STStimOut = dict(zip(range(self.nChannels), [SpikeList([], [])
             for i in range(self.nChannels)]))
        t_start = 0
        t_stop = 0
        if ch_events.get_nev() > 0:
            if normalize:
                ch_events = self.normalizeAER(ch_events)
                t_start = 0
            else:
                t_start = ch_events.get_all_tm().min() * 1e-3
            t_stop = ch_events.get_all_tm().max() * 1e-3

        if t_start == t_stop:
            t_stop += 1

        for ch in ch_events:
            stas = self[ch]
            #To ms
            ch_events.set_tm(ch, ch_events.get_tm(ch) / 1000)
            #eventsLogicalChannel[:,0]=eventsChannel[channelIdx][:,0]
            STStimOut[ch] = SpikeList(
                    ch_events.get_adtmev(ch),
                    ch_events.get_ad(ch),
                    t_start=t_start,
                    t_stop=t_stop
                    )

        return STStimOut

    def rawoutput_from_chevents(self, ch_events, func=None, normalize=True, filter_duplicates=False):
        """
        this function acts like generateST, but constructs a RawOutput object which delays the decoding until it is necessary.
        *ch_events* is a channelEvents object of type 'p' (Physical)

        Inputs:
            *func* - a dictionary of functions with channels as keys to decode the addresses. If omitted, all the channels are considered
            Outputs a RawOutput object
        """
        t_start = 0
        t_stop = 0

        ch_events = ch_events.copy()

        if func == None:
            # If no decoding functions are defined then use all of those
            # available in the channel addressing
            func_data = dict(zip(range(self.nChannels), [self[i]
                .addrPhysicalLogical for i in range(self.nChannels)]))
        else:
            func_data = func

        assert ch_events.atype == 'p', 'channel events must be of type physical (p)'

        if ch_events.get_nev() > 0:
            if normalize:
                ch_events = self.normalizeAER(ch_events)
                t_start = 0
            else:
                t_start = ch_events.get_first_tm() * 1e-3
            t_stop = ch_events.get_last_tm() * 1e-3

        if t_start == t_stop:
            t_stop += 1

        raw_data = {}

        for ch in ch_events:        
            ch_events.set_tm(ch, ch_events.get_tm(ch))
            raw_data[ch] = ch_events[ch]
        raw_out = RawOutput(
                raw_data,
                func_data,
                t_start=t_start,
                t_stop=t_stop,
                filter_duplicates=filter_duplicates
                )

        return raw_out

    def exportAER(self, spikeLists, filename=None, format='a', isi=True, sep='\t', addr_format='%u', time_format='%u', *args, **kwargs):
        '''
        spikeLists can be of the follwing type:
        - Monitors
        - SpikeList
        - list of SpikeLists of dimension nChannels
        - dictionary with channels as keys and SpikeLists as values
        - SpikeList is given, it will be interpreted as {0: spikeLists}.
        format specifies whether timestamps (format='t') or addresses (format='a') should be on the first column.
        *addr_format* and *time_format* format to be used by np.savetxt
        '''

        out = []
        assert format in ['t', 'a'], 'Format must be "a" or "t"'

        if hasattr(spikeLists, 'to_chstlist'):
            #Assuming it is of type Monitors
            spikeLists = spikeLists.to_chstlist()
        elif isinstance(spikeLists, list):
            assert len(spikeLists) == self.nChannels, "spikeLists must have dimension %d" % self.nChannels
            for i in range(len(spikeLists)):
                if not isinstance(spikeLists[i], SpikeList):
                    raise TypeError(
                        "Elements of spikeLists must be SpikeList objects!")

            spikeLists = dict(zip(range(len(spikeLists)), spikeLists))
        elif isinstance(spikeLists, SpikeList):
            spikeLists = {0: spikeLists}
        elif isinstance(spikeLists, dict):
            for i in spikeLists.iterkeys():
                if not isinstance(spikeLists[i], SpikeList):
                    raise TypeError(
                        "Values of spikeLists must be SpikeList objects!")
        else:
            raise RuntimeError(
                "spikeLists must be either a: SpikeList, list or dict object")

        ev = events(atype='Physical')

        #Translate logical addresses to physical using a mapping
        tic = time.time()
        for ch in spikeLists:
            if isinstance(spikeLists[ch], SpikeList):
                slrd = spikeLists[ch].raw_data()
                if len(slrd) > 0:
                    #addr_present=spikeLists[ch].id_list()
                    # mapping=dict(zip(addr_present,self.addrPhysicalConstruct(
                    # {ch:self[ch].addrLogicalExtract(addr_present)})))
                    # Not using convert because it is very slow ( iterates over
                    # all events )
                    tmp_mapped_SL = np.fliplr(slrd)
                    mapped_SL = np.zeros_like(tmp_mapped_SL, dtype='uint32')
                    mapped_SL[:, 1] = tmp_mapped_SL[:, 1] * 1000  # ms
                    mapped_SL[:, 0] = self[ch].addrLogicalPhysical(tmp_mapped_SL[:, 0]) + self.getValue(ch)
                    ev.add_adtmev(mapped_SL)
                    # ev.add_adtmev(mapSpikeListAddresses(spikeLists[ch],mappin
                    # g).convert(format='[id,time*1000]'))
                else:
                    print("Warning: Empty SpikeList encountered")
        tictoc = time.time() - tic
        if PERFORMANCE_DEBUG:
            print("Address encoding took {0} seconds".format(tictoc))

        #Multiplex
        tic = time.time()
        sortedIdx = np.argsort(ev.get_tm())
        tictoc = time.time() - tic
        if PERFORMANCE_DEBUG:
            print("Multiplexing took {0} seconds".format(tictoc))

        tic = time.time()
        #Create new sorted events object
        if len(sortedIdx) > 0:
            ev = events(ev.get_adtmev()[sortedIdx, :], atype='p')
            #exportAER
            if isi:
                ev.set_isi()

        else:
            ev = events(atype='p')

        #Choose desired output: no filename given, return events
        if filename is None:
            return ev
        #Otherwise write to file
        else:
            #exportAER_file
            if format is 'a':
                np.savetxt(filename, ev.get_adtmev(),
                    fmt=addr_format + sep + time_format)
            elif format is 't':
                np.savetxt(filename, ev.get_tmadev(),
                    fmt=time_format + sep + addr_format)
            return ev

    def buildAllHashTables(self, channels=None):
        if channels is None:
            channels = self.channels
        for i in channels:
            addrBuildHashTable(self[i])

STChannelAddressing = channelAddressing


def addrLogicalConstruct(stas, addr):
    """
    Constructs Logical addresses, *i.e.* addresses in float format useful for plotting because they keep the neuron - synapse order.

    **NOTE:** Logical addresses do not contain channel information

    *addr*: address (in a form `isValidAddress <pyst.isValidAddress>`_ can understand)
    """

    #Check and Parse Address
    nEntries, addr = isValidAddress(stas, addr)
    #Initialize logical address vector (double)
    addrLogical = np.zeros([nEntries], 'float')

    #Integer component
    IntCmp = np.zeros([nEntries])
    #Fractional component
    FracCmp = np.zeros([nEntries])
    IntBitsUsed = np.zeros([nEntries])
    FracBitsUsed = np.zeros([nEntries])

    for hrf_index, hrf in enumerate(stas.iter_hrfs()):  # hrf_index
        #Sanity check
        assert np.all((addr[hrf_index, :] & (2 ** hrf['bits'] - 1))
             == addr[hrf_index, :]), "Cropped significant bits"

        if hrf['type'] is 1:
            IntCmp = IntCmp + addr[hrf_index, :] * 2 ** IntBitsUsed
            IntBitsUsed = IntBitsUsed + hrf['bits']

        elif hrf['type'] is -1:
            # Considered using length instead: creates rational numbers such as
            # 0.33333333 which does not look nice
            FracBitsUsed = FracBitsUsed + hrf['bits']
            FracCmp = FracCmp + addr[hrf_index, :] * 2 ** (-FracBitsUsed)

    addrLogical = np.array(IntCmp + FracCmp, 'float')
    return addrLogical


def addrLogicalExtract(stas, addrLogical):
    """
    Description here
    """
    #TODO: Check and Parse Address
    addrLogical = np.array(addrLogical, 'float')
    nEntries = addrLogical.shape[0]
    FracBits = stas.nbits[-1]
    IntBits = 0

    #assert nEntries>0
    #Initialize logical address vector (double)
    addr = np.zeros([stas.nDims, nEntries], 'uint32')

    addr_frac, addr_int = np.modf(addrLogical)
    addr_int = addr_int.astype('int')
    addr_frac = (addr_frac * 2 ** stas.nbits[-1]).astype('int')

    for hrf_index, hrf in enumerate(stas.iter_hrfs()):  # hrf_index
        if hrf['type'] is 1:
            mask = (2 ** hrf['bits'] - 1)
            addr[hrf_index, :] = (addr_int >> IntBits) & mask
            IntBits += hrf['bits']
        elif hrf['type'] is -1:
            FracBits -= hrf['bits']
            addr[hrf_index, :] = (addr_frac >> FracBits) & (
                2 ** hrf['bits'] - 1)
            #decimal part

    return addr


def addrPhysicalConstruct(stas, addr):
    """
    Description here
    """

    #Check and Parse Address
    nEntries, addr = isValidAddress(stas, addr)

    #Initialize physical address vector (integer)
    addr = stas.addr_encoder.encode(addr)

    addrPhysical = [[]] * len(stas.field)
    for fieldIndex, field in enumerate(stas.iter_fields()):
        addrPhysical[fieldIndex] = field.construct(addr[fieldIndex])
    addrPhysical = np.sum(zip(*addrPhysical), 1)

    #Build dictionary for quick reference
    stas.addrPhysicalExtract(addrPhysical)

    return np.array(addrPhysical, 'uint32')


def isValidPhysicalAddress(stas, addrPhys):
    """
    Checks if the given list, int, numpy array is a valid physical address and outputs a numpy array
    """

    #Parse type, output numpy array
    if      isinstance(addrPhys, int):
        addrPhys = np.array([addrPhys], 'uint32')
    elif isinstance(addrPhys, np.int):
        addrPhys = np.array([int(addrPhys)], 'uint32')
    elif isinstance(addrPhys, np.int32):
        addrPhys = np.array([int(addrPhys)], 'uint32')
    elif isinstance(addrPhys, np.int64):
        addrPhys = np.array([int(addrPhys)], 'uint32')
    elif isinstance(addrPhys, list):
        addrPhys = np.array(addrPhys, 'uint32')
    elif isinstance(addrPhys, np.ndarray):
        addrPhys = addrPhys.astype('uint32')
# assert addrPhys.dtype==np.dtype('uint32'), 'Error: physical addresses must be
# an integer array or an integer'
    else:
        addrPhys = np.array(addrPhys, 'uint32')
        raise TypeError("Type was %s" % type(addrPhys))

    if len(addrPhys) is 0:
        return np.array([], 'uint32')
    else:
        return addrPhys


def addrPhysicalExtract(stas, addrPhys):
    """
    addrPhysicalExtract takes a physical address as argument and returns a numpy array containing the human readable list of addresses. First dimension is the field type, second dimension is the address.
    Uses hash tables for accelerating look-up.
    """
    return addrPhysicalExtractDecode(stas, addrPhys)


def addrPhysicalLogical(stas, addrPhys):
    """
    Direct translation from Physical Addresses to Logical addresses using a hash table
    """
    addrPhys=isValidPhysicalAddress(stas, addrPhys)
    #failedIndex=np.where(stas.addrExtractLogicalFast[addrPhys]==-1)[0]
    failedIndex = np.setdiff1d(addrPhys, stas.addrExtractLogicalFast.keys())
    #Don't even bother calling decode if everyone is in
    if len(failedIndex)>0:
       try:
           addrPhysicalLogicalDecode(stas, failedIndex) #Failed? decode the address
       except AssertionError as e:
           print('Error in addrPhysicalLogicalDecode.')
           print('No. of addresses in the packet : {0}'.format(len(addrPhys)))
           raise e


# NOTE: If you get a ValueError, you probably have checkLevel=0 and forgot to
# build the hash tables using addrBuildHashTable
    fastAddr = np.array(
        map(stas.addrExtractLogicalFast.get, addrPhys), 'float')

    # Loop should run over addresses only once if all addresse are in dict.
    # Small overhead
    return fastAddr.transpose()


def addrLogicalPhysical(stas, addrLogical, *args, **kwargs):
    """
    Direct translation from Logical Addresses to Physical addresses using a hash table
    """
    #Don't even bother calling decode if everyone is in
    try:
        # Failed? decode the address
        return addrLogicalPhysicalDecode(stas, addrLogical)  
    except AssertionError as e:
        print('Error in addrLogicalPhysicalDecode.')
        print('No. of addresses in the packet : {0}'.
            format(len(addrLogical)))
        raise e



def _buildGrid(inlist):
    nD = len(inlist)
    min_list = [0] * nD
    max_list = [None] * nD
    for i in range(nD):
        max_list[i] = len(inlist[i])
    strmgrid = 'np.mgrid['
    for fieldIdx in xrange(nD):
        strmgrid = strmgrid + '{0}:{1},'.format(
            min_list[fieldIdx], max_list[fieldIdx])
    strmgrid = strmgrid + ']'
    allpo = eval(strmgrid)
    tot = np.prod(allpo[0].shape)
    grid = np.zeros([tot, len(allpo)], dtype=type(min_list[0]))
    grid_values = np.zeros_like(grid)
    for j in range(len(allpo)):
        grid[:, j] = allpo[j].flatten()
    for i in range(nD):
        grid_values[:, i] = np.array(inlist[i])[grid[:, i]]
    return grid_values


def addrBuildHashTable(stas):
    """
    addrBuildHashTable(stas) constructs all possible physical addresses of the given address specification and stores it internally in its hash table (stas.addrExtractFast). This function is typically used for real-time monitoring

    *stas:* addrSpec object.

    *Hash tables are not used anymore*

    This function does not scale well.
    """
    nDim = []
    for hrf in stas.iter_hrfs():
        nDim.append(hrf['range'])
    pass

    #Construct grid
    addrLens = [None] * len(stas)

    for i in range(len(stas)):
        addrLens[i] = stas[i]['range']

    grid = _buildGrid(addrLens)

    stas.allpos = stas.addrLogicalConstruct(grid.transpose())
    addr_phys = addrLogicalPhysicalDecode(stas, stas.allpos)
    addrPhysicalLogicalDecode(stas, addr_phys)


def addrPhysicalExtractDecode(stas, addrPhys):
    """
    addrPhysicalExtract takes a physical address as argument and returns a list containing the addresses in human readable form.
    This is called only the first time the physical address is decoded.

    *stas:* addrSpec object.

    *addrPhys*: an integer or an iterable containing integers representing physical addresses
    """

    #initalize address list
    if not hasattr(addrPhys, '__len__'):
        addrPhys = [addrPhys]
    addrPhys = addrPhys & (2 ** stas.nBitsTotal - 1)

    addr = np.zeros([len(stas.field), len(addrPhys)], 'uint32')
    for fieldIndex, field in enumerate(stas.iter_fields()):
        addr[fieldIndex, :] = field.extract(addrPhys)

    addr = stas.addr_encoder.decode(addr)

    return addr


def addrPhysicalLogicalDecode(stas, addrPhys):
    '''
    Takes physical addresses and fills the physical to logical hash table
    '''
    addr = addrPhysicalExtractDecode(stas, addrPhys)
    addrLog = addrLogicalConstruct(stas, addr)
    #stas.addrExtractLogicalFast[addrPhys] = addrLog
    stas.addrExtractLogicalFast.update(zip(addrPhys, addrLog))
    return addrLog


def addrLogicalPhysicalDecode(stas, addrLog):
    '''
    Takes logical addresses and fills the logical to physical hash table
    '''
    addr = addrLogicalExtract(stas, addrLog)
    addrPhys = addrPhysicalConstruct(stas, addr)
    stas.addrExtractPhysicalFast.update(zip(addrLog, addrPhys))
    return addrPhys


def isValidAddress(stas, addrList):
    """
    This is an internal function which verifies that the given "human readable" address is consistent with the address specification. It also takes care of "filling in" the addresses: for example [range(15),2] is understood as [range(15),[2]*15].

    Raises AssertionError if addrList has incorrect data

    *stas*: an address specification object `addrSpec <pyST.addrSpec>`_
    *addrList*: a list of human readable addresses such as [range(15), 5]. It also accepts numpy arrays and transforms it accordingly
    """

        #    return len(addrList[0]),addrList
    if not isinstance(addrList, np.ndarray):
        assert getattr(addrList, '__iter__',
             False), "The addrList argument must be an iterable"
        # At the most minimum, check that the number of fields is consistent
        # with the length of the list
        assert len(addrList) == stas.nDims, "the number of fields is not consistent with the length of the list, are you using the correct address specification?"
        addrList = list(addrList)
    #make sure we have a 2 dimensional list and modify if necessary
        dimAddr = [0] * stas.nDims
        for hrf_index in stas.iter_hrf_index():
            if not hasattr(addrList[hrf_index], '__len__'):
                addrList[hrf_index] = [int(addrList[hrf_index])]
                #Make a list out of integer
            dimAddr[hrf_index] = len(addrList[hrf_index])

        #Remove all occurences of 1 in dimension list
        for i in xrange(dimAddr.count(1)):
            dimAddr.remove(1)

#        if not (len(np.unique(dimAddr)) in [0,1]):
# print("Dimensionalites of Address list is not consistent: Constructing grid")

        #Determine the number of entries
        if dimAddr != []:
            nEntries = np.prod(map(len, addrList))
        else:
            nEntries = 1

        #Initialize address list
        addrListFilled = [[-1] * nEntries for i in xrange(stas.nDims)]

        addrListFilled = _buildGrid(addrList).astype('uint32').transpose(
            )  # not good make buildGrid use same type as in the min/max lists

    else:
        addrListFilled = addrList.astype('uint32')
        nEntries = addrListFilled.shape[1]

    #Range Check
    for hrf_index, hrf in enumerate(stas.iter_hrfs()):
        assert np.all(addrListFilled[hrf_index, :]
            >= 0), "Address must be a postive integer."
        diffs = np.setdiff1d(
            np.unique(addrListFilled[hrf_index, :]), hrf['range'])
        if len(diffs) != 0:
            try:
                err_ind = np.where(addrListFilled[hrf_index] == diffs)[0]
            except Exception as e:
                print(np.where(addrListFilled[hrf_index] == diffs), diffs)
                raise AssertionError(None, None)
            wrongaddr = addrListFilled[:, err_ind]
            print("Address {3} is not in Range list ({0},{1}). Offending addresses on dimension {2}.".format(
                    np.max(hrf['range']),
                    np.max(hrf['range']),
                    hrf_index,
                    wrongaddr))
            raise AssertionError(wrongaddr, err_ind)
    return nEntries, addrListFilled


# NOTE: Scipy weave isn't ported to python 3
try:
    from scipy import weave
    from scipy.weave import converters
except:
    import weave
    from weave import converters
# ...


def _extract(x, a, r):
    '''
    Internal functions for applying the bit permutations defined in the chip file
    '''
    N = int(len(x))
    d = int(r.shape[0])
    r2 = 2 ** r
    a2 = 2 ** a
    y = np.zeros([N], 'int32')

    code = """
           for (int i=0; i<N; ++i) {
               for (int j=0; j<d; ++j) {
                    y[i]+=((x[i] & (a2[j]))>>a[j])*r2[j];
               }
           }
           """
    # compiler keyword only needed on windows with MSVC installed
    ext = weave.inline(code,
                       ['x', 'N', 'd', 'a', 'a2', 'r2', 'r', 'y'],
                       compiler='gcc')
    return y


def _construct(x, a, r):
    '''
    Internal functions for applying the bit permutations defined in the chip file
    '''
    N = int(len(x))
    d = int(r.shape[0])
    r2 = 2 ** r
    a2 = 2 ** a
    y = np.zeros([N], 'int32')

    code = """
           for (int i=0; i<N; ++i) {
               for (int j=0; j<d; ++j) {
                    y[i]+=((x[i] & (r2[j]))>>r[j])*a2[j];
               }
           }
           """
    # compiler keyword only needed on windows with MSVC installed
    ext = weave.inline(code,
                       ['x', 'N', 'd', 'a', 'a2', 'r2', 'r', 'y'],
                       compiler='gcc')
    return y


######################
######### STAS #######
class addressEncoder:
    def __init__(self, addr_conf, addr_str, addr_pinconf):
        addr_conf = _process_addrConf(addr_conf)
        addr_spec, nBits, nBitsTotal = _stas_parse_addrstr(addr_str)

        #Construct Extract, Construct
        self.fc, fc_field = _stas_create_construct(addr_conf, addr_pinconf)
        self.fe, fe_field = _stas_create_extract(addr_conf, addr_pinconf)

        nDims = len(addr_conf)

        self.fc_field_dict = {}
        self.fe_field_dict = {}

        id_list = extract_id_list(addr_pinconf)
        for i in range(len(addr_pinconf)):
            self.fc_field_dict[id_list[i]] = fc_field[i]

        id_list = extract_id_list(addr_conf)
        for i in range(len(addr_conf)):
            self.fe_field_dict[id_list[i]] = fe_field[i]

    def __getitem__(self, item):
        if item in self.fc_field_dict.keys():
            return self.fc_field_dict[item]
        elif item in self.fe_field_dict.keys():
            return self.fe_field_dict[item]
        else:
            raise AttributeError("There is no such field or pin %s" % item)

    def encode(self, addr):
        return np.vstack(self.fc(addr))

    def decode(self, addr):
        return np.vstack(self.fe(addr))

    def encode_field(self, field, addr):
        return np.array(self[field](*addr), 'uint32')

    def decode_field(self, field, addr):
        return np.array(self[field](*addr), 'uint32')


class layoutFieldEncoder:  # private
    """
    Internal struct used for addrSpec, performing the pin layout permutation (the information contained in addrSr)
    """

    def __init__(self, aspec, nWidth, position=0, pin=''):
        self.aspec = aspec
        self.aspec2 = 2 ** aspec
        self.nWidth = nWidth
        self.rWidth = np.arange(nWidth)
        self.r2Width = 2 ** self.rWidth
        self.position = position
        self.pin = pin
        #Check if transformation is an identity function. If so, just mask and shift
        if np.all(aspec == np.arange(aspec[0],aspec[-1]+1)):
            mask = np.sum([2**i for i in aspec]).astype('uint32')
            self.extract = lambda x: (x&mask)>>aspec[0]
            self.construct = lambda x: (x)<<aspec[0]
        else:
            self.extract = lambda x: _extract(x, self.aspec, self.rWidth)
            self.construct = lambda x: _construct(x, self.aspec, self.rWidth)


#TODO: Create fieldstruct for logical part and get rid of addrConf
#
#

#def generateST(id_list, events, normalize=True):
#    """
#    Extracts events from eventsChannel, using the address specification specified for the given channels, returns a list of SpikeList
#    Inputs:
#    *stas*: STas.addrSpec object
#    *events*: STas.events object
#    *normalize*: if True (default), the timestamp of the first event is set to zeros and all the timestamps of the subsequent events are shifted accordingly.
#    """
#    assert events.atype == 'l', 'channelEvents must be of type logical'
#
#    spike_list = SpikeList([], [])
#    t_start = 0
#    t_stop = 0
#
#    if normalize:
#        events.normalize_tm()
#
#    if events.get_nev() > 0:
#        t_start = events.get_tm().min() * 1e-3
#        t_stop = events.get_tm().max() * 1e-3
#
#    if t_start == t_stop:
#        t_stop += 1.
#
#    events.set_tm(events.get_tm() * 1e-3)
#    spike_list = SpikeList(
#            spikes=events.get_adtmev(),
#            id_list=id_list,
#            t_start=t_start,
#            t_stop=t_stop
#            )
#
#    return spike_list

TYPE_TO_NAME_DICT = { -1 : 'synapse dimension', 0 : 'other dimension', 1 : 'neuron dimension'}

def repr_addr_spec(addr_spec, nBitsTotal):
    '''
    Prints the address specification in a human readable format
    *Input*: addr_spec, a AddrSpec object
    '''
    base_string = '| ' + ''.join(['{'+str(i)+'} ' for i in range(nBitsTotal)[::-1]]) +' |'
    tmp_str = [ 'I ']*nBitsTotal
    for k,a in addr_spec.addrSpec.iteritems():
        a_rev = a[::-1] #reverse
        for i in a_rev:
            tmp_str[i] = '{0}{1}'.format(k,str(i))
    p2a_map = base_string.format(*tmp_str)
    base_string = ['']
    tmp_str = '{0} {1} = {2} \n'
    for n, a in enumerate(addr_spec.addrPinConf):
        field = addr_spec[n]
        base_string.append( tmp_str.format( TYPE_TO_NAME_DICT[field['type']], field['id'], field['f']) )
    a2n_map = ''.join(base_string) + 'I = ignore \n'
    return p2a_map + ' \n' + a2n_map + ' \n'



class addrSpec:
    """
    Address specification class
    # NOTE : Documentation is outdated!!

    *addrConf*: should be a list of dicts. Each dict provides the address specificatoin information for one field, and should contain the following entries:

    - 'f' : Function for translating the address of the field into its physical counterpart
    - 'id' : whose value is equal to a single, unique, lowercase character.
    - 'range' : an iterable containing the possible values the address field can take
    - 'type': which is either 1 (Neuron) or -1 (Synapse) or -2 (Connection parameter).

    *addrStr*: should contain a list of space-separated pins with a number following an uppercase character defining the position on the chip. For example: "X0 X1 X2 X3 X4 X5 X6 Y4 Y3 Y2 Y1 Y0"

    """
    def __init__(self, addrStr='', addrConf=list(), addrPinConf='', id="NoName", nhml = False):
        ##script for parsing the adddrStr
        self.id = id
        self.addrPinConf = addrPinConf
        self.addrStr = addrStr
        self.addrConf = addrConf
        #Useful to have the following functions as members:
        self.__class__.addrLogicalConstruct = addrLogicalConstruct
        self.__class__.addrLogicalExtract = addrLogicalExtract
        self.__class__.addrPhysicalConstruct = addrPhysicalConstruct
        self.__class__.addrPhysicalExtract = addrPhysicalExtract
        self.__class__.addrPhysicalLogical = addrPhysicalLogical
        self.__class__.addrLogicalPhysical = addrLogicalPhysical
        self.__class__.repr_addr_spec = repr_addr_spec
        #Update all parameters
        if not nhml:
            self.update()

    def update(self):
        '''
        (Re)Generates the object based on the addrStr, addrConf and addrPinConf
        '''

        self.addrConf = _process_addrConf(self.addrConf)
        self.addrDict = _process_addrDict(self.addrConf)
        self.addrSpec, self.nBits, self.nBitsTotal = _stas_parse_addrstr(
            self.addrStr)
        self.addr_encoder = addressEncoder(
            self.addrConf, self.addrStr, self.addrPinConf)
        self.nDims = len(self.addrConf)
        self.field, self.nFields = _stas_create_fields(
            self.nBits, self.addrSpec, self.addrConf, self.addrPinConf)
        self.nbits = _stas_compute_nbits(self.addrConf)
        #self.addrExtractLogicalFast = np.empty(
        #    [2 ** np.sum(self.nbits.values())], 'float')
        # NOTE: The above code is modified to accomodate for blank bits in the
        # address space and hence use of nBitsTotal might be more accurate.
        self.addrExtractLogicalFast = dict()
        self.addrExtractPhysicalFast = dict()
#        try:
        # Building addresses on the fly.
        # Cuses memory errors otherwise for large AER spaces
        #addrBuildHashTable(self)
#        except Exception as e:
#            warnings.warn('Could not BuildHashTable: {0}'.format(e))

    def __len__(self):
        return self.nDims

    def __repr__(self):
        input_format = "[ "
        for v in self.iter_hrfs():
            input_format += str(v['id']) + ", "
        input_format += "]"

        return "AddrSpec: " + self.id + " " + input_format

    def iter_hrf_index(self):
        '''
        Iterates over the indexes of addrConf
        '''
        for k in xrange(self.nDims):
            yield k

    def iter_hrfs(self):
        '''
        Iterates over the list addrConf
        '''
        for v in self.addrConf:
            yield v

    def iter_fields(self):
        for field in self.field:
            yield field

    def iter_fields_by_pin(self, pin_order=None):
        if pin_order == None:
            pin_order = np.sort(self.nBits.keys())
        for p in pin_order:
            for field in self.field:
                if field.pin == p:
                    yield p, field

    def __getitem__(self, field):
        return self.addrConf[field]

    def __getNHML__(self):
        '''
        Returns xml representation of this object
        '''
        doc = etree.Element('addressSpecification')
        # addrConf
        for conf in self.addrConf:
            # Dimensions
            dimdoc = etree.SubElement(doc, 'dim')
            # Id
            dimdoc.attrib['id'] = conf['id']
            # Type
            if conf['type'] == 1:
                dimdoc.attrib['type'] = 'soma'
            elif conf['type'] == -1:
                dimdoc.attrib['type'] = 'synapse'
            elif conf['type'] == -2:
                dimdoc.attrib['type'] = 'connection'
            else:
                # NOTE: In future if there are more types this is where they
                # are
                # to be added!
                dimdoc.attrib['type'] = 'unknown'
            # Range
            rangedoc = etree.SubElement(dimdoc, 'range')
            rangedoc.text = str(conf['range'])
            # Description
            descdoc = etree.SubElement(dimdoc, 'description')
            # Decoder
            decoderdoc = etree.SubElement(dimdoc, 'decoder')
            decoderdoc.text = conf['f']
        # Pins
        for pins in self.addrPinConf:
            pindoc = etree.SubElement(doc, 'pin')
            pindoc.attrib['id'] = pins['id']
            decoderdoc = etree.SubElement(pindoc, 'decoder')
            decoderdoc.text = pins['f']
        # Pinlayout
        pinlayout = etree.SubElement(doc, 'pinlayout')
        pinlayout.text = self.addrStr
        return doc

    def __parseNHML__(self, doc):
        '''
        Parses an lxml element tree or a file name with xml content to
        initialize the object.
        '''
        if isinstance(doc, str):
            # parse the file
            doc = etree.parse(doc).getroot()
        else:
            # assuming doc is an lxml Element object.
            assert doc.tag == 'addressSpecification'

        self.addrConf = []  # NOTE: should probably be done in init function
        self.addrPinConf = []  # NOTE: should also be done in init function
        self.addrStr = ''  # NOTE: should also be done in init function
        # addrConf
        for elm in doc:
            if elm.tag == 'dim':
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
                        if chld.text:
                            dim['range'] = eval(chld.text)
                        else:
                            dim['range'] = None
                    elif chld.tag == 'description':
                        # Description
                        dim['description'] = chld.text
                    elif chld.tag == 'decoder':
                        # Decoder
                        dim['f'] = chld.text
                    else:
                        pass
                self.addrConf.append(dim)
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
            elif elm.tag == 'pinlayout':
                # Pinlayout
                self.addrStr = elm.text
            else:
                pass
        # Update the object
        self.update()
        return


#################################################################
#The following functions are used in AddrSpec
#
def _process_addrConf(addrConf):
    #Find number of required bits from range data
    for hrf in addrConf:
        bits_req = 1
        while max(hrf['range']) >> bits_req > 0:
            bits_req += 1
        hrf['bits'] = bits_req
        hrf['default_value'] = 0  # TODO: Add it to chipfile
    return addrConf


def _process_addrDict(addrConf):
    d = dict()
    #Find number of required bits from range data
    for i, hrf in enumerate(addrConf):
        d[hrf['id']] = i
    return d


def _stas_create_fields(nBits, addrSpec, addrConf, addrPinConf):
    """
    Creating fields for physical
    """
    nFields = len(nBits)
    field = [None for i in range(nFields)]
    pos = 0

    for nextfield in addrPinConf:
        aspec = addrSpec[nextfield['id']]
        currentfield = layoutFieldEncoder(
                        aspec=np.array(aspec, 'uint'),
                        nWidth=nBits[nextfield['id']],
                        position=pos,
                        pin=nextfield['id'],
                        )
        field[pos] = currentfield
        pos += 1

    return field, nFields


def _stas_compute_nbits(addrConf):
    nbits = {-1: 0, 1: 0, 0: 0, -2:0}
    for a in addrConf:
        nbits[a['type']] += a['bits']
    return nbits


def extract_id_list(addr_conf):
    id_list = []
    for i in addr_conf:
        id_list.append(i['id'])
    return id_list


def _stas_create_extract(addr_conf, addr_pinconf):
    """
    Parsing the Extraction functions
    """

    arg_str = ''.join([i + ', ' for i in extract_id_list(addr_pinconf)])
    fe_field = [eval('lambda ' + arg_str + ':' + s['f']) for s in addr_conf]
    fe = eval('lambda ' + arg_str + ':' + "[" + ''.join([s['f'] +
        ', ' for s in addr_conf]) + ']')
    nDims = len(addr_pinconf)  # Trick to find out d
    return lambda x: fe(*x), fe_field


def _stas_create_construct(addr_conf, addr_pinconf):
    """
    Parsing the Construction functions
    """

    arg_str = ''.join([i + ', ' for i in extract_id_list(addr_conf)])
    fc_field = [eval('lambda ' + arg_str + ':' + v['f']) for v in addr_pinconf]
    fc = eval('lambda ' + arg_str + ':' + "[" + ''.join([s['f'] +
        ', ' for s in addr_pinconf]) + ']')
    return lambda x: fc(*x), fc_field


def _stas_parse_addrstr(addrStr):

    assert isinstance(addrStr, str), "addrStr must be a string!"
    addrSpec = {}
    nBits = {}

    addrStrsplit = addrStr.split()
    addrStrsplit.reverse()

    for i in addrStrsplit:
        if i[0] in addrSpec:
            if len(addrSpec[i[0]]) < (int(i[1:]) + 1):
                addrSpec[i[0]].extend([None] * (int(i[1:]) + 1 -
                    len(addrSpec[i[0]])))
        else:
            addrSpec[i[0]] = [None] * (int(i[1:]) + 1)
            nBits[i[0]] = 0

        addrSpec[i[0]][int(i[1:])] = addrStrsplit.index(i)
    nBitsTotal = 0
    for k, v in addrSpec.items():
        nBitsTotal += len(v)
        if k == 'I': #I is the ignore bit
            addrSpec.pop(k)
            nBits.pop(k)
        else:
            nBits[k] = len(v)

    #nBitsTotal = np.sum(nBits.values())
    return addrSpec, nBits, nBitsTotal


def load_stas_from_csv(CSVfile):
    """
    This is a reduced version of pyNCS.Chip._readCSV, which reads only the addressing information
    Inputs:
    *CSVfile* - A csv file describing the chip (a "chipfile").
    """
    if CSVfile == None:
        raise IOError('filename must be given')

    def indexfrom(list, str):
        '''Return list index whcih contains specified string'''
        for i in range(len(list)):
            if str in list[i].lower():
                return i

    minwords = ['channel', 'signal', 'bias']
    keywords = ['channel', 'block', 'signal', 'bias', 'range',
        'description', 'fet', 'pad', 'value']
    validtable = 0

    csv = None
    with open(CSVfile, 'r') as CSV:
        csv = CSV.readlines()
    cb = ''
    #flags
    validtable = False
    stasStim_f = False
    stasMon_f = False
    addrConf = []
    addrPinConf = []
    addrStr = ''

    aerIn, aerOut = None, None
    for line in csv:
        line = line.replace('\'', '')  # remove single and double quotes
        line = line.replace('"', '')
        if 'chipclass' in line.lower():  # Chipclass
            chipclass = line.strip().split('\t')[1]
        elif keywords[0] in line.lower():  # Table header
            fields = line.strip().split('\t')
            for i in minwords:
                if not indexfrom(fields, i):
                    raise Exception('Bad file format: missing required column : %s, requires channel, signal and bias' % i)
            validtable = True
        elif line.startswith('pinlayout') and (stasStim_f or stasMon_f):  # addStr specification
            addrStr = line.strip().split('\t')[1]
        elif line.startswith('aerIn'):
            stasStim_f = True
            addrConf = []
            fc_str = ''
            fe_str = ''
        elif line.startswith('aerOut'):
            stasMon_f = True
            addrConf = []
            addrPinConf = []
            fc_str = ''
            fe_str = ''
        elif line.startswith('id') and (stasStim_f or stasMon_f):
            x = line.strip().split('\t')
            d = {}
            d['id'] = x[1].strip()
            d['range'] = eval(x[3])
            d['range_str'] = x[3]
            d['type'] = eval(x[5])
            d['f'] = x[7]
            addrConf.append(d)

        elif line.startswith('pinid') and (stasStim_f or stasMon_f):
            x = line.strip().split('\t')
            d = {}
            d['id'] = x[1].strip()
            d['f'] = x[3]
            addrPinConf.append(d)
        elif line.startswith('\t') or line.startswith('\n'):
            validtable = False
            if stasStim_f:
                aerIn = addrSpec(addrStr=addrStr, addrConf=addrConf,
                     addrPinConf=addrPinConf, id=chipclass + "In")
                stasStim_f = False
            elif stasMon_f:
                aerOut = addrSpec(addrStr=addrStr, addrConf=addrConf,
                     addrPinConf=addrPinConf, id=chipclass + "Out")
                stasMon_f = False

    #Determining the dimensions from the address specifications
    return aerIn, aerOut

def load_stas_from_nhml(doc):
    '''
    load_stas_from_nhml(NHMLfile)
    This function returns addrSpec objects aerIn and aerOut by parsing the NHML
    file.
    '''
    if isinstance(doc, str):
        # parse the file
        doc = etree.parse(doc).getroot()
    else:
        # assuming doc is an lxml Element object
        assert doc.tag == 'chip'
    chipclass = doc.get('chipclass')
    aerIn, aerOut = None, None
    for elm in doc:
        if elm.tag == 'addressSpecification':
            if elm.get('type') == 'aerIn':
                aerIn = addrSpec(id=chipclass + 'In', nhml = True)
                aerIn.__parseNHML__(elm)
            elif elm.get('type') == 'aerOut':
                aerOut = addrSpec(id=chipclass + 'Out', nhml = True)
                aerOut.__parseNHML__(elm)
            else:
                pass
    return aerIn, aerOut
    
