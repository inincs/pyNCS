#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Fabio Stefanini
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
from __future__ import absolute_import
import numpy as np
import copy
import itertools as it
from pickle import dump, load

from .group import AddrGroup


def _buildGrid(inlist):
    '''
    Builds a multidimensional grid from input lists.
    (tip: operate on indexes, not values)
    Ex:
    >>> _buildGrid([[0,1],[0,2,3]])
    >>> array([[0, 0],
               [0, 2],
               [0, 3],
               [1, 0],
               [1, 2],
               [1, 3]])
    Sadique: Isn't this is the same as
    np.array(list(it.product(*inlist))) ??
    '''
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


# TODO: Population.soma[0] should give the (Human) address
class Population(object):
    """
    Population is a set of neurons and corresponding synapses. Population can have parameters
    (efficacy of excitatory input synapses, ...).
    This is on top of synapses and is intended to be used by the user to create neural networks.
    """
    # TODO: extend to multiple chips
    # TODO: proper pickling

    def __init__(self, name = '', description = '',
                 setup=None, chipid=None, neurontype=None):
        """
        Init a population by name and description. Population is empty.
        Name and description are used in the graph representation for
        connectivity.
        - name: string
        - description: string
        - setup: NeuroSetup to init the population
        - chipid: string id of the chip in the setup, e.g. 'ifslwta'
        - neurontype: neurontype string to init the population, e.g. 'pixel'
        """
        self.name = name
        self.description = description
        self.soma = AddrGroup('Empty group.')
        self.synapses = dict()
        self.setup = None
        self.neuronblock = None
        if setup is not None and chipid is not None and neurontype is not None:
            self.__populate_init__(setup, chipid, neurontype)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value
        if hasattr(self, 'soma'):
            self.soma.name = self._name + ' ' + 'soma'

        if hasattr(self, 'neuronblock'):
            for s in self.neuronblock.synapses.keys():
                self.synapses[self.neuronblock.synapses[s].
                    id].name = self._name + ' ' + str(s)

    def __copy__(self):
        '''
        Return an exact copy of the population.
        '''
        p = Population(self.name, self.description)
        p.soma = copy.deepcopy(self.soma)
        p.synapses = copy.deepcopy(self.synapses)
        p.setup = self.setup
        p.neuronblock = self.neuronblock
        return p

    def __getitem__(self, i):
        """
        x.__getitem__(i) ==> x[i]
        """
        p = Population(self.name, self.description)
        p.setup = self.setup
        p.neuronblock = self.neuronblock

        if isinstance(i, slice):
            p.soma = self.soma[i]
        else:
            if hasattr(i, '__len__'):
                p.soma = self.soma[i]
            else:
                p.soma = self.soma[i:i + 1]
        p.__populate_synapses__()
        return p

    def __getslice__(self, i, j):
        """
        x.__getslice__(i, j) ==> x[i:j]
        """
        p = Population(self.name, self.description)
        p.setup = self.setup
        p.neuronblock = self.neuronblock
        p.soma = self.soma[i:j]
        p.__populate_synapses__()
        return p

    def __len__(self):
        """
        len(x) <==> len(x.soma)

        Return the number of neurons in the population.
        """
        return len(self.soma)

    def clear(self):
        """
        Clear the population back to its original state.
        """
        self.soma = AddrGroup('Empty group')
        self.synapses = dict()

#    def __populateByExplicitAddr__(self, chipid, addr):
#        """
#        This function is useful if you know the addresses of neurons and
#        synapses. Needs consistence between chipid and addr.
#        WARNING: you are supposed to use higher level functions, so use this at
#        your own risk!
#        """
#        raise NotImplementedError

    def isinit(self):
        """
        Return True if population is initiated with setup, chipid and
        neurontype, e.g. not populated.
        """
        if self.setup is None or self.neuronblock is None:
            return False
        else:
            return True

    def union(self, population):
        """
        Add the given population's addresses to the existing one. If the
        address is already there, it doesn't add.
        WARNING: Synapses might be incorrect
        """
        if not population.isinit():
            # if the other population is still not initiated init with the same
            # parameters of self population but don't populate
            population.init(self.setup, self.soma.chipid, self.neuronblock.id)
        if population.neuronblock.neurochip.id\
           != self.neuronblock.neurochip.id:
                raise Exception(
                    'Union between different chips is not yet implemented.')
        # TODO: actually, addr should be constructed with the appropriate
        # dimensions

        if len(self.soma.addr) != 0 and len(population.soma.addr) != 0:
            self.soma.addr = np.hstack([self.soma.addr,
                                           population.soma.addr])
        else:  # self.soma is empty
            self.soma.addr = population.soma.addr.copy()

        self.soma.repopulate(self.setup)

        # At this point the order of the addresses is dictated by the order by
        # which they have been added
        # Also, synapses order may also be broken, consider adding
        # __populate_synapses__() here (too slow)
        self.soma.sort()

        return

    def add(self, addresses):
        """
        Adds a neuron (with all its synapses) to the population. Population has to
        be populated already. Address has to be of the appropriate format for
        the chip on which the population has been allocated.
        Arguments are:
            - addresses: neuron address in human format (e.g. [10, 2] for neuron
              [10,2] in a 2D chip.
        """
        for a in addresses:
            if not a in self.neuronblock.soma.addresses:
                raise Exception("At least one address is not present in the population neuronblock.")

        # add soma address(es)
        self.soma.add(self.setup, addresses)
        # add synapses
        for k in self.synapses.keys():
            ch = self.synapses[k].channel
            soma_addr = addresses
            syn_addr = self.__soma2syn__(soma_addr, synapses=[k])
            self.synapses[k].add( self.setup, syn_addr)
        #self.__unique_addresses__()

    def remove(self, address):
        """
        Removes a neuron (with all its synapses) from the population.
        Arguments are:
            - address: neuron address in human format (e.g. [10, 2] for neuron
              [10,2] in a 2D chip.
        """
        # self.soma.remove([address])...
        raise NotImplementedError

    def __getstate__(self):
        """
        Implement pickling functionalities when dumping.
        """
        d = dict(self.__dict__)
        d['neurontype'] = self.neuronblock.id
        del d['neuronblock']  # don't need to dump neuronblock
        return d

    def __setstate__(self, dict):
        """
        Implement pickling functionalities when loading.
        """
        # TODO: can we do this with __getinitargs__ instead?
        self.__dict__ = copy.copy(dict)
        chip = dict['setup'].chips[dict['soma'].chipid]
        neuronblock = chip.neuron[dict['neurontype']]
        self.neuronblock = neuronblock
        del self.neurontype

    def __soma2syn__(self, addresses, synapses=None):
        """
        Given the neurons addresses, returns all the addresses of its synapses.
        Useful for populating by giving only the soma addresses.
        If synapses is set, returns only those synapses (e.g., 'excitatory0').
        NOTE: The returned array preserves ordering of neurons (np.repeat)
        """

        somaaddr = addresses
        dtp = AddrGroup._get_dtype(AddrGroup('',''), self.setup, self.soma.chipid, 'in')
        syn_dim_names = None
        synaddrs_all = None
        for s in synapses:
            syn_block = self.neuronblock.synapses[s]
            if syn_dim_names == None:
                syn_dim_names = syn_block.dims.keys()
                syndtp = np.dtype([(str(k), 'uint32') for k in syn_dim_names])
                # NOTE: the str type cast is because numpy dtype doesn't support
                # unicode strings
            synaddrs = [syn_block.dims[fld] for fld in syn_dim_names]
            synaddrs = list(it.product(*synaddrs))
            if synaddrs_all == None:
                synaddrs_all = np.array(synaddrs)
            else:
                synaddrs_all = np.concatenate((synaddrs_all, synaddrs))
        # Syn addresses only
        synaddrs_all = synaddrs_all.astype('uint32').view(syndtp).T[0]
        # Syn addresses tiled for all soma addresses
        synaddrs_full = np.tile(synaddrs_all, len(addresses)) 
        # Soma addresses repeated for all synapse addresses
        soma_addr_all = addresses.repeat(len(synaddrs_all))
        # Full addresses space
        syn_addr_full = np.zeros(soma_addr_all.shape, dtype=dtp)
        for fld in soma_addr_all.dtype.names:
            syn_addr_full[fld] = soma_addr_all[fld]
        for fld in synaddrs_all.dtype.names:
            syn_addr_full[fld] = synaddrs_full[fld]
        return syn_addr_full

    def __populate_init__(self, setup, chipid, neurontype):
        """
        Basic operations common to every populate method.
        """
        self.setup = setup

        chip = setup.chips[chipid]

        # check whether neuron is available in the chip
        try:
            neuronblock = chip.neuron[neurontype]
        except KeyError:            
            raise KeyError('ERROR: %s: No such neurontype in current setup.' % neurontype)

        self.neuronblock = neuronblock  # neuronblock contains translations
                                        # for somas AND synapses biases

        self.soma.__populate__(setup, chipid, 'out')
        self.soma.name = self.name + ' ' + 'soma'

    def init(self, setup, chipid, neurontype):
        """
        self.init ==> self.__populate_init__
        """
        return self.__populate_init__(setup, chipid, neurontype)

    def __populate_synapses__(self):
        """
        Populate all the synapses of the population with the corresponding
        addresses for the neurons in the population.
        """
        #self.soma.sort() #sadique: is this necessary ??
        for s in self.neuronblock.synapses.keys():
            syn_id = self.neuronblock.synapses[s].id
            self.synapses[syn_id] = S = AddrGroup(self.
                neuronblock.synapses[s].id)
            self.synapses[syn_id].name = self.name + ' ' + str(s)

            ##Populate empty first to set channel and ch_addr of S
            ##Consider doing this in AddrGroup
            #S.populate_line(self.setup,
            #        self.soma.chipid,
            #        grouptype='in',
            #        addresses=[])
            addresses=self.__soma2syn__(self.soma.addr, [s])
            S.populate_line(self.setup,
                    self.soma.chipid,
                    grouptype='in',
                    addresses=addresses,)

    def populate_all(self, setup, chipid, neurontype):
        """
        Populate all the neurons in the given chip.
        """
        self.__populate_init__(setup, chipid, neurontype)

        # filter addresses from the ones available in neuronblock
        addresses = self.neuronblock.soma.addresses
        self.soma.populate_line(setup, chipid, grouptype='out',
                                addresses=addresses)

        self.__populate_synapses__()

    def populate_by_number(self, setup, chipid, neurontype, n, offset=0):
        """
        Takes the given number of addresses from the neuronblock available
        addresses. It takes the first n addresses if offset is not set. Arguments
        are:
            - setup: a NeuroSetup
            - chipid: id of the chip as expressed in setup.xml
            - neurontype: id of neurons as expressed in chipfile.xml (e.g.
            'excitatory')
            - n: the number of neuron to allocate
            - offset: imposes not to take the first addresses
        """

        self.__populate_init__(setup, chipid, neurontype)

        # filter addresses from the ones available in neuronblock
        addresses = self.neuronblock.soma.addresses[offset:n + offset]
        if len(addresses) != n:
            raise Exception("Not enough neurons of this type.")
        self.soma.populate_line(setup, chipid, grouptype='out',
                                addresses=addresses)

        self.__populate_synapses__()

    def populate_sparse(self, setup, chipid, neurontype, p=.3):
        """
        Populate picking random addresses from the neuronblock of all possible
        addresses with probability p. Arguments are:
            - setup: a NeuroSetup
            - chipid: id of the chip as expressed in setup.xml
            - neurontype: id of neurons as expressed in chipfile.xml (e.g.
            'excitatory')
            - p: probability of picking neurons in [0, 1)
        """

        self.__populate_init__(setup, chipid, neurontype)

        # filter addresses from the ones available in neuronblock
        a = np.array(self.neuronblock.soma.addresses)
        addresses = a[np.random.random(len(a)) < p]

        self.soma.populate_line(setup, chipid, grouptype='out',
                                addresses=addresses)

        self.__populate_synapses__()

    def populate_by_topology(self, setup, chipid, neurontype, topology='rectangle', topology_kwargs={'p1': [0, 0], 'p2': [63, 0]}):
        """
        Takes the given number of addresses by respecting the chips natural topology (i.e. 0 to n%X+n/Y). It takes the first n addresses if offset is not set. Arguments
        are:
            - setup: a NeuroSetup
            - neurontype: id of neurons as expressed in chipfile.xml (e.g.
            'excitatory')
            - n: the number of neuron to allocate
        """

        self.__populate_init__(setup, chipid, neurontype)

        # all the addresses for the given topology
        S = AddrGroup(self.neuronblock.soma.id)
        getattr(S, 'populate_' + topology)(setup, chipid,
             grouptype='out', **topology_kwargs)
        addresses = _set_intersection(S.addr, self.neuronblock.soma.addresses)
        self.soma.populate_line(
            setup, chipid, grouptype='out', addresses=addresses)
        self.__populate_synapses__()

#    def populate_by_dimension(self, setup, chipid, neurontype, filt_list = [[-1],[-1],[-1]]:
#        """
#        """
#
#        self.__populate_init__(setup, chipid, neurontype)
#
#        def f(x):
#            if all([xx in filt_list[i] for i,xx in enumerate(x)])
#            
#            
#

    
    def populate_by_addr_list(self, setup, chipid, neurontype, id_list=[]):
        """
        Assigns the population a set of human readable addresses explicityly
        delared by the user.
        Arguments are:
            - setup: a NeuroSetup
            - chipid: id of the chip as expressed in setup.xml
            - neurontype: id of neurons as expressed in chipfile.xml (e.g.
            'excitatory')
            - id_list: the list of ids for neurons to allocate (human readable addresses).  

        *Example*: populate 3 neurons on a 3d grid x,y,z with addresses
        [1,0,3], [1,2,3] and [2,2,2]

        .. code-block:python

        >> populate_by_idlist(setup, chipid, neurontype, 
                              id_list = [[1,0,3],[1,2,3],[2,2,2]])
        """
        self.__populate_init__(setup, chipid, neurontype)
        try:
            self.soma.populate_line(
                setup, chipid, grouptype='out', addresses=id_list)
        except Exception as e:
            print("Chip {0} does not contain one or more of the specified addresses.".format(chipid))
            raise e
        self.__populate_synapses__()

    
    def populate_by_id(self, setup, chipid, neurontype, id_list=[], axes=[]):
        """
        Look documentation for populate_by_filter
        Function depricated, use populate_by_addr_list or populate_by_idfilter in
        future.
        """
        self.populate_by_id_filter(setup, chipid, neurontype,
                                   id_filter=id_list, axes=axes)
        
    def populate_by_id_filter(self, setup, chipid, neurontype, id_filter=[], axes=[]):
        """
        Filters addresses from the available neuronblock addresses. 
        Arguments are:
            - setup: a NeuroSetup
            - chipid: id of the chip as expressed in setup.xml
            - neurontype: id of neurons as expressed in chipfile.xml (e.g.
            'excitatory')
            - axes: list of axes by which to filter the addresses
            - id_filter: the list of ids (per dimension) for neurons to allocate (human readable addresses).  

        *Example*: populate all neurons on a 3d grid x,y,z such that 24<= x < 26, y=5:

        .. code-block:python

        >> populate_by_id_filter(setup, chipid, neurontype, id_filter =
        [[24,25], [5]], axes = [0,1])

        where x is defined to be on axis = 0, and y to be on axis = 1
        """

        self.__populate_init__(setup, chipid, neurontype)
        #Backward compatibility
        if not hasattr(axes, '__iter__'):
            axes = [axes]
        
        # Ensuring all elements in the filter list are lists
        for i, l in enumerate(id_filter):
            if not hasattr(l, '__iter__'):
                id_filter[i] = [l]

        #filter addresses fast , but uses slightly more memory
        mask = np.zeros_like(self.neuronblock.soma.addresses, dtype='bool')
        for i, a in enumerate(axes):
            mask[:,a]=np.in1d(self.neuronblock.soma.addresses[:,a], id_filter[i])
        addresses = self.neuronblock.soma.addresses[np.prod(mask, axis=1, dtype='bool')]
        try:
            self.soma.populate_line(
                setup, chipid, grouptype='out', addresses=addresses)
        except:
            raise Exception(("Chip {0} contains no neurons of given" +
                            "id_list on axes {1}.").format(chipid, axes))

        self.__populate_synapses__()


def _set_intersection(A, B):
    """
    returns elements that are both in A and B
    Does NOT recover original order!!! Instead sorts according to the same order as logical addresses
    """
    A = set(tuple(i) for i in A)
    B = set(tuple(i) for i in B)
    AinterB = np.array(list(A.intersection(B)), 'int')
    #Importing from unordered list: must sort back to initial order
    # zip(*a[::-1]) because we want to sort with respect to rightmost column
    # (this is how pyST builds it)
    #AinterB=AinterB[np.lexsort(zip(*AinterB[:,::-1]))]
    return _sort_by_logical(AinterB)


def _sort_by_logical(S):
    """
    Sort in-place according to the logical addesses order
    """
    return S[np.lexsort(zip(*S[:, :]))]


def _sort_by_antilogical(S):
    """
    Sort in-place according to the logical addesses order
    """
    return S[np.lexsort(zip(*S[:, ::-1]))]


def _flatten(l, ltypes=(list, tuple, np.ndarray)):
    '''
    The function flattens lists of lists or tuples or ints
    Taken from web
    http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
    '''
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            # The commented code is supposed to eliminate blank ararys
            # But has a bug and also eliminates zeros
            #if not l[i] and l[i]!=0:
            #    l.pop(i)
            #    i -= 1
            #    break
            #else:
            l[i:i + 1] = l[i]
        i += 1
    return ltype(l)
