#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Fabio Stefanini
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
from __future__ import absolute_import
import itertools
from warnings import warn
from .pyST import *

class AddrGroupBase(object):
    def __init__(self, name, description=None):
        if description:
            self.description = description
        else:
            self.description = 'group'
        self.name = name  # a short name
        self._channel = None
        self.chipid = ''
        self.addr = None
        self._paddr = None
        self._laddr = None
        self._channel = None
        self.grouptype = ''
        self.dtype = None
    
    def __deepcopy__(self, memo):
        return self.__copy__()

    def __copy__(self):
        '''
        Return an exact copy of the group.
        '''
        g = AddrGroup(self.name, self.description)
        g._channel = self._channel
        g.chipid = self.chipid
        g.addr = self.addr.copy()
        g.setup = self.setup
        g.addrspec = self.addrspec
        g.grouptype = self.grouptype
        g.dtype = self.dtype
        g.repopulate()
        return g

    @property
    def channel(self):
        return self._channel
    
    @channel.setter
    def channel(self, ch):
        self._channel = ch 
    
    @property
    def laddr(self):
        return self._laddr
    
    @property
    def paddr(self):
        return self._paddr   
    
    @laddr.setter
    def laddr(self, laddr):
        self._laddr = laddr
    
    @paddr.setter
    def paddr(self, paddr):
        self._paddr = paddr    
 
     
    def __len__(self):
        """
        x.__len__() ==> len(x)
        """
        return len(self.addr)

    def __eq__(self, other):
        '''
        Compares two address groups
        '''
        return (self.paddr == other.paddr).all()

    def __getitem__(self, i):
        """
        x.__getitem__(i) ==> x[i]
        """
        g = self.__copy__()
        if isinstance(i, slice):
            g.addr = self.addr[i].copy()
            if self._paddr != None:
                g._paddr = self._paddr[i]
            if self._laddr != None:
                g._laddr = self._laddr[i]
        else:
            if hasattr(i, '__len__'):
                g.addr = np.array(self.addr[i])
                if self._paddr != None:
                    g._paddr = np.array(self._paddr[i], dtype='uint32')
                if self._laddr != None:
                    g._laddr = np.array(self._laddr[i], dtype='float')
            else:
                g.addr = np.array([self.addr[i]], dtype=self.dtype)
                if self._paddr != None:
                    g._paddr = np.array([self._paddr[i]], dtype='uint32')
                if self._laddr != None:
                    g._laddr = np.array([self._laddr[i]], dtype='float')
        return g

    def __getslice__(self, i, j):
        """
        x.__getslice__(i, j) ==> x[i:j]
        """
        g = self.__copy__()
        g.addr = self.addr[i:j]
        if self._paddr != None:
            g._paddr = self._paddr[i:j]
        if self._laddr != None:
            g._laddr = self._laddr[i:j]
        return g


    def sort(self, order=None):
        """
        Sort all the addresses with the given order.
        If the order is None it orders by the last column of the address (most probably the column of synapses).
        Otherwise order must be a list of the fields in self.addr. 
        This is very useful for sorting 2D+ addresses.        
        """
        if order is None:
            idx = np.argsort(self.addr)
        else:
            idx = np.lexsort(map(self.addr.__getitem__, order))
        # Just to ensure the addresses are generated.
        self.laddr; self.paddr
        # Sort the addresses accordingly.
        self.addr = self.addr[idx]
        self._laddr = self._laddr[idx]
        self._paddr = self._paddr[idx]

    def add(self, addresses):
        """
        Add a list of addresses to a group.
        """
        if self.grouptype is '':
            raise Exception("Group has never been populated before.")
        
        if len(addresses) == 0:
            return

        # Check for data type
        if isinstance(addresses, np.ndarray):
            if np.issubdtype(addresses.dtype, self.dtype):
                addr = addresses
            else:
                # WARNING: This data type conversion can introduce errors for
                # complex data types.
                addr = addresses.astype('uint32').view(self.dtype)
                addr = addr.T[0]
        else:
            addr = np.array(addresses, dtype='uint32').view(self.dtype)
            addr = addr.T[0]

        if len(self.addr) > 0:
            self.addr = np.concatenate([self.addr, addr])
        else:
            self.addr = addr


    def is_empty(self):
        return len(self) == 0

    # TODO: convert spiketrains methods to spiketrains(fashion='poisson') style
class AddrGroup(AddrGroupBase):
    #TODO: merge function
    """
    A AddrGroup is a group of chip addresses. They can be of grouptype 'in' or
    of grouptype 'out'. A description is needed for populating the group.

    Ex.:
        setup = pyNCS.Setup('setuptype.xml')
        setup.load('setup')
        setup.apply()

        # output from the retina
        group1 = AddrGroup('retina output')
        group1.populate_rectangle(
            setup, 'retina', 'out', [0,0,0], [64,64,0], n = 28*28 )

        # input from the ifslwta_0 chip
        group2 = AddrGroup('A1 cortex learning input')
        group2.populate_rectangle(setup, 'ifslwta_0', 'in', [0,4], [28,32])

        # input from the sequencer
        group3 = AddrGroup('A1 cortex excitatory input')
        group3.populate_rectangle(setup, 'sequencer', 'in', [0,2], [28,2])

        # input from the sequencer
        group4 = AddrGroup('A1 cortex excitatory inhibitory')
        group4.populate_rectangle(setup, 'ifslwta_0', 'in', [0,2], [28,2])

        group5 = AddrGroup('56 trains from sequencer')
        group4.populate_by_number(setup, 'seq', 'out', 56)

        # you want to treat group3 and group4 as a whole
        group4.join_with(group3)

        mapping = []
        mapping = group1.connect_one2one(group2)
             # the connect functions will maybe go in some Mapping module
        # ... and so on and so forth ...

        setMappings(mapping)

    """



    @property
    def ch_addr(self):
        if self.grouptype == 'in':
            return self.setup.seq
        elif self.grouptype == 'out':
            return self.setup.mon
        elif self.grouptype == '':
            warn('Grouptype has not been set')
            return None
        elif self.grouptype == '':
            raise ValueError('Grouptype should be None, "in" or "out", not {0}'.
                format(self.grouptype))

    @property
    def channel(self):
        if self._channel == None:
            # if it is None generate the necesary addresses
            if self.grouptype == 'in':
                self._channel =\
                self.setup.aerSlot[self.setup.chipslots[self.chipid]]['seqIn'][0]
            elif self.grouptype == 'out':
                self._channel =\
                self.setup.aerSlot[self.setup.chipslots[self.chipid]]['monOut'][0]

        return self._channel



    def repopulate(self, setup = None):
        """
        Repopulates laddr and paddr with respect to addr
        """
        if setup is not None:
            self.setup = setup.chaddrspecs
        if self.grouptype is '':
            raise Exception("Group has never been populated before.")
        self._channel = None
        self._laddr = None
        self._paddr = None
        self.addr = self.addr.copy()


#    def remove(self, address):
#        """
#        Remove list of addresses to a group.
#        """
#
#        if self.grouptype is '':
#            raise Exception("Group has never been populated before.")
#
#        raise NotImplementedError

    def __populate__(self, setup, chipid, grouptype, addresses=[]):

        self.chipid = chipid
        self.grouptype = grouptype
        self.setup = setup.chaddrspecs

        if self.grouptype is 'in':
            self.addrspec = setup.seq[self.channel]
        elif self.grouptype is 'out':
            self.addrspec = setup.mon[self.channel]
        else:
            self.addrspec = None

        self.dtype = self._get_dtype(setup, chipid, grouptype)
        self.addr = np.array([], dtype=self.dtype)
        self._paddr = None  # np.array([],dtype='uint32')
        self._laddr = None  # np.array([],dtype='float')
        self.add(addresses)
    
    def _get_dtype(self, setup, chipid, grouptype):
        '''
        Returns the data type of AddrGroup.addr variable, ie. the
        format of human readable addresses.
        '''
        if grouptype == 'in':
            flds = setup.get_chip_aerin(chipid).addrDict
        elif grouptype == 'out':
            flds = setup.get_chip_aerout(chipid).addrDict
        else:
            raise ValueError('Grouptype should be None, "in" or "out", not {0}'.
                format(grouptype))
        dtp = [None for i in range(len(flds))]
        for k,v in flds.iteritems():
            dtp[v] = (k, 'uint32')
        return np.dtype(dtp)


    def populate_line(self, setup, chipid, grouptype, addresses):
        self.__populate__(setup, chipid, grouptype, addresses)

    def populate_rectangle(self, setup, chipid, grouptype, p1, p2, z=None, n=None):
        """
        Populate with all the neurons within the rectangle's given coordinates. If n is given,
        populate with the given number of neurons instead of all neurons in the rectangle.
        z defines the 3d coordinate, e.g. retina polarity or synapse in the 2d chip.
        """

        if n:
            shaperatio = (1. * p2[1] - p1[1]) * 1. / (1. * p2[0] - p1[0])
            nx, ny = int(
                np.sqrt(1. * n / shaperatio)), int(np.sqrt(1. * shaperatio * n))
            addrx = np.linspace(p1[0], p2[0], nx)
            addry = np.linspace(p1[1], p2[1], ny)
        else:
            addrx = np.arange(p1[0], p2[0] + 1)
            addry = np.arange(p1[1], p2[1] + 1)

        if z == None:
            addresses = np.column_stack([addrx.repeat(len(
                addry)), np.array(addry.tolist() * len(addrx))])
        else:
            addresses = np.column_stack([addrx.repeat(len(addry)), np.array(addry.
                tolist() * len(addrx)), np.repeat(z, len(addrx) * len(addry))])

        self.__populate__(setup, chipid, grouptype, addresses)

    def populate_cuboid(self, setup, chipid, grouptype, p1, p2):
        '''
        Populate addresses within volume defined by the two opposite vertices
        in the N-dimensional space.
        setup: setup object
        chipid: ID of the chip to be populated on
        grouptype: IN or OUT
        p1: Lower Vertex (included) (vector of size n) defining the starting vertex
        p2: Upper Vertex (excluded) (vector of size n) defining the end vertex
        '''
        try:
            assert( len(p1)==len(p2) )
        except AssertionError as e:
            raise Exception("Dimensions of vertices do not match")
        try:
            assert( ((np.array(p2) - np.array(p1)) <= 0).sum() == 0 )
        except AssertionError as e:
            warn('Zero volume for given vertices')
        
        edges = []
        for i in range(len(p1)):
            edges.append(range(p1[i], p2[i]))
        addresses = itertools.product(*edges)
        self.__populate__(setup, chipid, grouptype, list(addresses))

#    def populate_by_number(self, setup, chipid, grouptype, n, dims):
#        """
#        Populate with the given number of neuron, doesn't matter the shape.
#        Needs the number of dimensions
#        WARNING: This should go in upper layer of abstraction.
#        TODO: check if this function is unused (Reason: contains reference errors)
#        """
#        max = np.prod([len(i) for i in dims])
#        if n > max:
#            raise Exception(
#                'Don\'t have enough space on this chip (max is %d).' % max)
#
#        lendims = [len(i) for i in dims]
#        nd = 0  # number of needed dimensions
#        div = 1
#        while True:
#            div = int(n / np.prod(lendims[0:nd + 1]))
#            if div > 0:
#                print "Dimension %d is not enough" % nd
#                nd += 1
#            else:
#                break
#
#        # unfold dimensions here
#        self.__unfold_dims__(dims)
#
#        dim = len(xdims)
#
#        n_filled_lines = int(n / dim)
#        filled_lines = np.arange(n_filled_lines)
#
#        n_return_line = n % dim
#        return_line = np.arange(n_return_line)
#
#        addrx = np.array(xdims)
#
#        addresses = np.column_stack([addrx.repeat(n_filled_lines),
#             np.array(filled_lines.tolist() * dim)])
#        addresses_ret = np.column_stack(
#            [return_line, np.repeat(filled_lines, n_return_line)])
#
#        addresses = np.concatenate([addresses, addresses_ret])
#
#        self.__populate__(setup, chipid, grouptype, addresses)


    def spiketrains(self, spikes=[], t_start=0, t_stop=None, dims=None):
        '''
        Constructs a stimulus with the desired SpikeList.
        spikes : list of tuples (id, time)
        id -> logical addresses of the group ie. AddrGroup.laddr
        All the arguments are the same as that of a pyST.SpikeList
        '''
        stStim = SpikeList(spikes=spikes, id_list=self.laddr,
             t_start=t_start, t_stop=t_stop, dims=dims)
        return {self.channel: stStim}

    def spiketrains_single(self, spike_time, t_after=100):
        """
        Creates a spiketrain with a single event with the given spike_time for all the addresses in the group.
        Inputs:
        *spike_time*: a float or a vector of floats indicating the spike times. If a vector is given, each entry is used once for each address in the logical addresses of the group (in the same order as self.laddr)
        *t_after* ms after the spike
        """
        s = np.empty([len(self), 2])
        s[:, 1] = spike_time
        s[:, 0] = self.laddr
        t_stop = max(s[:, 1]) + t_after
        spiketrain = self.spiketrains(s, t_start=0, t_stop=t_stop)
        return spiketrain

    def spiketrains_regular_gaussian(self, rate, offset=0., scale=5., t_start=0., duration=1000., channel=None):
        '''
        Returns a SpikeTrain whose spikes are regularly spaced, but jittered
        according to a Gaussian distribution around the spiking period, with
        the given rate (Hz) and stopping time t_stop (milliseconds).

        Note: t_start is always 0.0, thus all realizations are as if they
        spiked at t=0.0, though this spike is not included in the SpikeList.

        Inputs:
            rate    - the rate of the discharge (in Hz). Can be an iterable
            t_start - the beginning of the SpikeTrain (in ms)
            phase   - Use an offset vector to add an offset shift to the spike trains. Can be an iterable
            scale   - width of the Gaussian distribution placed at the regular
                      spike times, according to which the spike will be drawn (in ms)
            duration- The duration of the SpikeTrain (in ms)
            channel - Channel argument to enforce a particular channel
        '''
        stStim = SpikeList([], id_list=[])

        if not hasattr(rate, '__len__'):
            rate = rate * np.ones_like(self.laddr)
        if not hasattr(offset, '__len__'):
            offset = offset * np.ones_like(self.laddr)
        elif len(self.laddr) != len(rate):
            raise RuntimeError('Rate vector must be of the same length as the number of neurons in population: {0}'.format(len(self.addr)))

        if self.is_empty():
            print("AddrGroup is empty!")
            return stStim

        for i, id in enumerate(self.laddr):
            stStim[id] = STCreate.regular_gaussian_generator(rate[i],
                                                    phase=offset[i],
                                                    scale=scale,
                                                    t_start=t_start,
                                                    t_stop=t_start + duration)

        if channel is None:
            channel = self.channel

        return {channel: stStim}

    def spiketrains_regular(self, rate, offset=0., t_start=0., duration=1000.,
                            jitter=False, channel=None):
        """
        Create a regular spiketrain for each address on the group.
        You can use channel argument to enforce a particular channel.
        Use an offset vector to add an offset shift to the trains (must be a
        vector).

        Inputs:
            rate    - the rate of the discharge (in Hz). Can be an iterable
            t_start - the beginning of the SpikeTrain (in ms)
            offet   - Offset the spiketrain by this number (in ms.). Can be an iterable.
            jitter  - whether the spiketrain should be jittered by an amount numpy.random.rand()/rate
            duration- The duration of the SpikeTrain (in ms)
            channel - Channel argument to enforce a particular channel
        """
        stStim = SpikeList([], id_list=[])

        if not hasattr(rate, '__len__'):
            rate = rate * np.ones_like(self.laddr)
        if not hasattr(offset, '__len__'):
            offset = offset * np.ones_like(self.laddr)
        elif len(self.laddr) != len(rate):
            raise RuntimeError('Rate vector must be of the same length as the number of neurons in population: {0}'.format(len(self.addr)))

        if self.is_empty():
            print("AddrGroup is empty!")
            return stStim

        for i, id in enumerate(self.laddr):
            stStim[id] = STCreate.regular_generator(rate[i],
                                                    phase=offset[i],
                                                    jitter=jitter,
                                                    t_start=t_start,
                                                    t_stop=t_start + duration)

        if channel is None:
            channel = self.channel

        return {channel: stStim}

    def spiketrains_poisson(self, rate, t_start=0., duration=1000.,
                            channel=None):
        """
        Create a poisson spiketrain for each address on the group.
        You can use channel argument to enforce a particular channel.
        """
        stStim = SpikeList([], id_list=[])

        if not hasattr(rate, '__len__'):
            rate = rate * np.ones_like(self.laddr)
        elif len(self.laddr) != len(rate):
            raise RuntimeError('Rate vector must be of the same length as the number of neurons in population: {0}'.format(len(self.addr)))

        if self.is_empty():
            print("AddrGroup is empty!")
            return stStim

        for i, id in enumerate(self.laddr):
            stStim[id] = STCreate.poisson_generator(
                rate[i], t_start=t_start, t_stop=t_start + duration)

        if channel is None:
            channel = self.channel

        return {channel: stStim}

    def spiketrains_inh_poisson(self, rate, t, channel=None, **kwargs):
        """
        Create inhomogeneous poisson spiketrains. Rate is a vector of rates.
        The first dimension of rates corresponds to the neurons
        The second dimension corresponds to the time bin. The length of this must be the same as t.
        In addition, the rates *must* end with 0 to mark the end of the last bin.
        keyword arguments kwargs are passed to the spiketrain generator.
        See also pyST.STCreate.inh_poisson_generator
        """
        return self.spiketrains_inh_generator(rate, t, channel, base_generator=STCreate.poisson_generator, **kwargs)

    def spiketrains_inh_generator(self, rate, t, channel=None, base_generator=None, **kwargs):
        """
        Create inhomogeneous spiketrains.
        The process is defined by the base_generator function. See pyST.STCreate for available functions.
        Rate is a vector of rates.
        The first dimension of rates corresponds to the neurons
        The second dimension corresponds to the time bin. The length of this must be the same as t.
        In addition, the rates *must* end with 0 to mark the end of the last bin.
        keyword arguments kwargs are passed to the spiketrain generator.
        See also pyST.STCreate.inh_poisson_generator
        Kewword arguments are passed to the spike generator. Default is regular.
        """
        if base_generator == None:
            base_generator = STCreate.regular_generator
        stStim = SpikeList([], id_list=[])

        if len(self.addr) != rate.shape[0]:
            raise RuntimeError('Rate vector must be of the same length as the number of neurons in population: {0}'.format(self.laddr))

        assert rate.shape[1] == len(
            t), "time vector must be compatible with time axis of rate"

        if self.is_empty():
            print("AddrGroup is empty!")
            return stStim
        for i, id in enumerate(self.laddr):
            stStim[id] = STCreate.inh_poisson_generator(rate[i],
                 t, t_stop=t[-1], base_generator=base_generator, **kwargs)

        if channel is None:
            channel = self.channel

        return {channel: stStim}

    def spiketrains_1D_bump(self, pos, width, ampl, t_start=0., duration=                            1000., channel=None):
        rate = self.__spiketrains_periodic_bump_vector(pos, width, ampl)

        return self.spiketrains_poisson(rate, t_start, duration,
                                        channel=channel)



    def __spiketrains_periodic_bump_vector(self, pos, width, ampl):
        from scipy.stats import norm
        N = len(self.addr) / 2
        rate = norm.pdf(np.arange(len(self.addr)), scale=width, loc=N)
        rate = np.roll(rate, pos - N)
        rate = rate / rate.max() * ampl
        return rate

    def getLaddr(self):
        '''
        Regenreates Locigal addresses from address list.
        '''
        if self._laddr != None:
            return self._laddr
        if len(self.addr) > 0:
            self._laddr = self.ch_addr[self.channel].addrLogicalConstruct(
                self.addr.view('uint32').reshape((-1, len(self.dtype))).T)
        else:
            self._laddr = np.array([], dtype='float')
        return self._laddr

    def setLaddr(self, ad):
        self._laddr = ad
        warn('Not advisable to modify AddrGroup.laddr. Change AddrGroup.addr and repopulate() instead')

    def getPaddr(self):
        '''
        Generates Physical addresses from the address list. To regenerate use self.repopulate
        '''
        if self._paddr != None:
            return self._paddr
        if len(self.addr) > 0:
            self._paddr = self.ch_addr.addrPhysicalConstruct(
                {self.channel: self.addr.view('uint32').reshape((-1,
                                                                 len(self.dtype))).T})
        else:
            self._paddr = np.array([], dtype='uint32')
        return self._paddr

    def setPaddr(self, ad):
        self._paddr = ad
        warn('Not advisable to modify AddrGroup.paddr. Change AddrGroup.addr and repopulate() instead')

    # Properties of AddrGroup
    laddr = property(getLaddr, setLaddr, doc='Logical addresses')
    paddr = property(getPaddr, setPaddr, doc='Physical addresses')
