#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Fabio Stefanini
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
import os
import sys
import subprocess
import numpy as np
from warnings import warn


class Mapping(object):
    """
    A class representing the mapping between groups of chip addresses.
    """

    def __init__(self, name, description=None):
        """
        """
        # TODO: proper pickling
        # TODO: restore graph from mappings

        if description:
            self.description = description
        else:
            self.description = 'This mapping describes the AER connectivity\
                                between addresses.'
        self.name = "\"" + 'AER connectivity' + "\""
        self.mapping = []

    def __getstate__(self):
        """
        """
        d = dict(self.__dict__)
        return d

    def __setstate__(self, dict):
        """
        """
        self.__dict__ = dict

    def __graph_from_mapping__(self):
        raise NotImplementedError

    def __len__(self):
        return len(self.mapping)

    def __toint__(self):
        """
        Convert all the elements to int because some numpy functions can return
        longint.#
        """
        self.mapping = [[int(x) for x in v] for v in self.mapping]

    def __instance_from_matrix_random(self, M):
        try:
            return np.random.binomial(1, M.tolist()).astype('bool')
        except:
            return np.random.binomial(1,M).astype('bool')

    def import_from_connections(self, connections_list):
        """
        Reads the list of connections and updates its mapping table.
        Connections must be pyNCS.Connection instancies.
        """
        self.clear()
        self.mapping = np.concatenate([c.mapping.mapping for c in
                                       connections_list]).tolist()

    def complete(self, connlist):
        '''
        Completes the fields if there are any blank fields in the table
        '''
        return connlist

    def merge(self, pyncs_mapping):
        """
        Merge the existing mapping with a given one.
        """
        # TODO: check for duplicates
        self.mapping = self.complete(self.mapping)
        if len(self.mapping) > 0:
            # check for shape
            my_cols = np.shape(self.mapping)[1]
            if len(pyncs_mapping.mapping) > 0:
                # check for shape
                their_cols = np.shape(pyncs_mapping.mapping)[1]
                # check for shape compatibility
                if my_cols > their_cols:
                    # Do something
                    self.mapping.extend(pyncs_mapping.mapping)
                    self.mapping = self.complete(self.mapping)
                elif my_cols < their_cols:
                    warn('Connection cannot be merged. Ignoring probability')
                    # NOTE: Assuming the only missing dimension is probability
                    self.mapping = np.concatenate([self.mapping,
                        np.array(pyncs_mapping.mapping)[:,0:my_cols]]).tolist()
                else:
                    self.mapping = np.concatenate([self.mapping,
                                        pyncs_mapping.mapping]).tolist()
        else:
            self.mapping = pyncs_mapping.mapping

    def connect(self, groupsrc, groupdst, expand=True, fashion='one2one', fashion_kwargs={}, connection_kwargs={}, check=True):
        """
        Wrap the connect call to all type of different connectivity functions.
        
        *fashion*: Connection fashion, arguments to specific fashion can be passed through keyword arguments *fashion_kwargs*.
                   Default fashion is 'one2one'.
        
        When check is True synapses cannot be output addresses and somas cannot be input addresses.

        *Example:*
        >>  mymapping.connect(src, dst, "all2all", expand = false)
        which is equivalent to:
        >>  mymapping.__connect_all2all__(src, dst, expand=False}
        """
        if not expand:
            self.clear()


        if check and not self.is_connect_possible(groupsrc, groupdst):
            return []

        #__connect__ function return a list of indexes of addresses
        connlist =  getattr(self, '__connect_' + fashion + '__', )(
                range(len(groupsrc)), 
                range(len(groupdst)),
                                                   **fashion_kwargs)

        #... and the following function creates a list of physical addresses
        self.mapping += self.connlist_to_mapping(groupsrc, groupdst, connlist, connection_kwargs)

        return self.mapping
        

    def connlist_to_mapping(self, groupsrc, groupdst, connlist, connect_kwargs):
        '''
        Creates a list of connections of physical addresses given the list of indexes in connlist.
        '''
        if len(connlist)>0:
            self.groupsrc = groupsrc.__copy__()[np.array(connlist)[:,0]]
            self.groupdst = groupdst.__copy__()[np.array(connlist)[:,1]]

            if connect_kwargs:
                for k,v in connect_kwargs.iteritems():
                    field_idx = self.groupdst.addrspec.addrDict[k]
                    #Slow, consider vectorizing
                    if not hasattr(v, '__len__'):
                        v = [v]*len(self.groupdst.addr)
                    for i, a in enumerate(self.groupdst.addr):
                        a[field_idx] = v[i]
                        self.groupdst.addr[i] = a
                    self.groupdst.repopulate()

            return np.column_stack([self.groupsrc.paddr,self.groupdst.paddr]).tolist()
        else:
            return []

        

    def __connect_all2all__(self, groupsrc, groupdst):
        """
        Connect groups in an all to all fashion.
        """
        return self.__connect_random_all2all__(groupsrc, groupdst, p=1.0)

    def __connect_one2one__(self, groupsrc, groupdst, p=1.0):
        """
        Connects in a one to one fashion, from the first of the source to the
        last of the source.  If the sizes are different it just raises a
        warning!
        """

        if len(groupsrc) != len(groupdst):
            print("WARNING: source and destination have different sizes")
            if len(groupdst) > len(groupsrc):
                groupdst = groupdst[:len(groupsrc)]
            else:
                groupsrc = groupsrc[:len(groupdst)]
        connect_dist = np.eye(len(groupsrc))*p
        return self.__connect_by_probability_matrix__(groupsrc, groupdst, connect_dist)

    def __connect_random_all2all__(self, groupsrc, groupdst, p=0.25):
        """
        Connects in an all to all fashion with probability p for every
        connection, from the first of the source to the last of the source.
        """

        connect_dist = np.ones([len(groupsrc), len(groupdst)]) * p

        return self.__connect_by_probability_matrix__(groupsrc, groupdst, connect_dist)

    def __connect_by_boolean_matrix__(self, groupsrc, groupdst, connection):
        '''
        groupsrc: source group
        groupdst: destination group
        connection: matrix of connections
        '''


        pairs_all = np.array(np.meshgrid(groupsrc, groupdst)).transpose()
        connection = np.array(connection, 'bool')
        pairs_selected = pairs_all[connection]

        c1 = pairs_selected[:, 0]
        c2 = pairs_selected[:, 1]

        return np.column_stack([c1, c2]).tolist()

    def __connect_by_probability_matrix__(self, groupsrc, groupdst, M, return_inst = False):

        connect_inst = self.__instance_from_matrix_random(M)
        mapping = self.__connect_by_boolean_matrix__(groupsrc, groupdst, connect_inst)

        if return_inst:
            return connect_inst
        else:
            return mapping

    def __resize_matrix_resample(self, groupsrc, groupdst, M):
        from scipy.signal import resample
        x_resampled = resample(M, len(groupsrc), window='blk')
        xy_resampled = resample(
            x_resampled.transpose(), len(groupdst), window='blk').transpose()
        return np.minimum(np.maximum(xy_resampled, 0), 1)

    def __connect_by_arbitrary_matrix__(self, groupsrc, groupdst, M, resize_method='resample'):

        M = getattr(self, '_Mapping__resize_matrix_' + resize_method)(
            groupsrc, groupdst, M)
        connect_inst = self.__instance_from_matrix_random(M)

        return self.__connect_by_boolean_matrix__(groupsrc, groupdst, connect_inst)

    def __connect_shuffle_all2all__(self, groupsrc, groupdst, p=0.25):
        """
        Connects in an all to all fashion by shuffling p*N/N 1's and 1-p*N/N 0's for every source, from the first of the source to the last of the source.
        Has the advantage that the number of connections is kept fixed, which can reduce inhomogeneity in the outputs.
        p can be an iterable, in which case the probability of each source to connect to the target i is equal to p[i].
        i.e. creates a matrix M such that P(M_{src,tgt}==1) = p[tgt] with src=1, ..., Nsrc and tgt=1, ... ,Ntgt
        """

        if not hasattr(p, '__len__'):
            p = p * np.ones(len(groupdst))

        N = [p[i] * len(groupsrc) for i in range(len(groupdst))]
        connect_inst = np.zeros([len(groupsrc), len(groupdst)]).astype('bool')

        # It is important to know that p is intepreted as the probability to
        # connect to the target, not that the source connects to any target
        #i.e. M such that P(M_{src,tgt}==1) = p[tgt]

        #Create N[i] for each source i...
        for i in range(len(groupdst)):
            #Separate integer and factional parts
            frac, integ = np.modf(N[i])
            connect_inst[:integ, i] = True
            # Deal with fractional part by rolling a dice and checking not to
            # overrun
            if np.random.rand() < frac and integ + 1 < len(groupsrc):
                connect_inst[integ + 1, i] = True

        #... and shuffle them
        for i in range(connect_inst.shape[1]):
            np.random.shuffle(connect_inst[:, i])

        #TODO: go through binary matrix

        return self.__connect_by_boolean_matrix__(groupsrc, groupdst, connect_inst)


    def clear(self):
        """
        Clear mapping list.
        """
        self.mapping = []


    def save(self, filename):
        """
        Save the mapping into a file.
        """
        np.savetxt(filename, self.mapping)


    def load(self, filename, verbose=False):
        """
        Loads the mapping from a file.
        """
        self.clear()
        self.mapping = np.loadtxt(filename)        


    def is_connect_possible(self, groupsrc, groupdst):
        if groupsrc.grouptype != 'out':
            warn("The source group is not of type 'out'")
        if groupdst.grouptype != 'in':
            warn("The target group is not of type 'in'")
        return True


    def prepare(self):
        '''
        Run any functions before applying the mapping table
        '''
        self.__toint__()
        
    def write(self, *args, **kwargs):
        warn("Mappings.write has no effect. Use NeuroSetupBase.run or NeuroSetupBase.prepare instead")

    def map_events(self, events):
        '''
        Map events in software.
        Input:
        *paddr*: a pyST.events object, containing physical addresses
        Output:
        evetns with mapped Physical addresses
        '''
        map_dict = self.mapping_dict()
        return events.filter_by_mapping(map_dict)

    def mapping_dict(self):
        '''
        Returns a dictionary of the mappings
        '''
        #sort list
        from collections import defaultdict
        mapping_dict = defaultdict(list)
        func = lambda srctgt: mapping_dict[srctgt[0]].append(srctgt[1])
        map(func, self.mapping)
        return mapping_dict


class PMapping(Mapping):
    """
    A class representing the mapping between groups of chip addresses, supporting connection probabilities
    """
    max_p = 127


    def __instance_from_matrix_random(self, M):
        return M

    def __connect_by_probability_matrix__(self, groupsrc, groupdst, M, return_inst=False):
        '''
        M : connection matrix, where each entry is the transmission probability
        '''

        #Create pairs of source - destinations and reshape as appropiately
        pairs = np.array(np.meshgrid(groupsrc, groupdst)).transpose().reshape(-1, 2)

        p = np.column_stack([pairs, M.flatten()])

        # Transform p. Somehow this is not the best place to do this since it
        # assumes already the discretization due to the mapper firmware
        p[:, 2] = np.array(p[:, 2] * self.max_p, 'uint32')

        #Remove all zero probabilities
        zp = p[:, 2] != 0
        p_nonzero = p[zp, :].astype('uint32').tolist()

        return p_nonzero


    def connlist_to_mapping(self, groupsrc, groupdst, connlist, connect_kwargs={}):
        '''
        Creates a list of connections of physical addresses given the list of indexes in connlist.
        '''
        if len(connlist)>0:
            self.groupsrc = groupsrc.__copy__()[np.array(connlist)[:,0]]
            self.groupdst = groupdst.__copy__()[np.array(connlist)[:,1]]

            connlist = self.complete(connlist)

            if connect_kwargs:
                for k,v in connect_kwargs.iteritems():
                    field_idx = self.groupdst.addrspec.addrDict[k]
                    #Slow, consider vectorizing
                    if not hasattr(v, '__len__'):
                        v = [v]*len(self.groupdst.addr)
                    for i, a in enumerate(self.groupdst.addr):
                        a[field_idx] = v[i]
                        self.groupdst.addr[i] = a
                    self.groupdst.repopulate()

            return np.column_stack([self.groupsrc.paddr,self.groupdst.paddr,np.array(connlist)[:,2]]).tolist()
        else:
            return []

    def complete(self, connlist):
        '''
        For all connections without a probability, assume the probability is one
        '''
        for c in connlist:
            if len(c) == 2:
                c.append(self.max_p)
        return connlist


    def prepare(self):
        super(PMapping,self).prepare()
