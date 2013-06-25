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
import pydot
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
        self.graph = pydot.Graph(self.name)

    def __getstate__(self):
        """
        """
        d = dict(self.__dict__)
        del d['graph']
        return d

    def __setstate__(self, dict):
        """
        """
        self.__dict__ = dict
        self.graph = pydot.Graph(self.name)

    def __graph_from_mapping__(self):
        pass

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

    def __build_from_pmatrix_shuffle(self, M):
        raise NotImplementedError

    def import_from_connections(self, connections_list):
        """
        Reads the list of connections and updates its mapping table.
        Connections must be pyNCS.Connection instancies.
        """
        self.clear()
        self.mapping = np.concatenate([c.mapping.mapping for c in
                                       connections_list]).tolist()

    def complete(self):
        '''
        Completes the fields if there are any blank fields in the table
        '''
        pass

    def merge(self, pyncs_mapping):
        """
        Merge the existing mapping with a given one.
        """
        # TODO: check for duplicates
        self.complete()
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
                    self.complete()
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

    def connect(self, groupsrc, groupdst, expand=True, fashion='one2one', fashion_kwargs={},
                check=True):
        """
        Wrap the connect call to all type of different connectivity functions.
        Arguments to specific fashion can be passed through keyword arguments.
        Default fashion is 'one2one'. When check is True synapses cannot be
        output addresses and somas cannot be input addresses.
        Example:
            mymapping.connect(src, dst, "all2all", {'expand':False}
        which is equivalent to:
            mymapping.__connect_all2all__(src, dst, expand=False}
        """
        if check and not self.is_connect_possible(groupsrc, groupdst):
            return []

        return getattr(self, '__connect_' + fashion + '__', )(groupsrc, groupdst,
                                                   expand=expand,
                                                   **fashion_kwargs)

    def __connect_all2all__(self, groupsrc, groupdst,
                            expand=True, hide=False):
        """
        Connect groups in an all to all fashion.
        """
        c1 = np.concatenate([np.repeat(i, len(groupdst.paddr)) for i in
                             groupsrc.paddr])
        c2 = list(groupdst.paddr) * len(groupsrc.addr)
        if expand:
            self.mapping.extend(np.column_stack([c1, c2]).tolist())
        else:
            self.mapping = np.column_stack([c1, c2]).tolist()

        if not hide:
            self.add_edge(groupsrc, groupdst, arrowhead='crow', dir='both')
        
        return self.mapping

    def __connect_one2one__(self, groupsrc, groupdst, p=1.0,
                            expand=True, hide=False):
        """
        Connects in a one to one fashion, from the first of the source to the
        last of the source.  If the sizes are different it just raises a
        warning!
        """

        if len(groupsrc) != len(groupdst):
            print "WARNING: source and destination have different sizes"
            if len(groupdst) > len(groupsrc):
                groupdst = groupdst[:len(groupsrc)]
            else:
                groupsrc = groupsrc[:len(groupdst)]
        connect_dist = np.eye(len(groupsrc))*p
        return self.__connect_by_probability_matrix__(groupsrc, groupdst,
                                                      connect_dist, expand=True, hide=False)

    def __connect_random_all2all__(self, groupsrc, groupdst, p=0.25,
                                   expand=True, hide=False):
        """
        Connects in an all to all fashion with probability p for every
        connection, from the first of the source to the last of the source.
        """

        connect_dist = np.ones([len(groupsrc), len(groupdst)]) * p

        return self.__connect_by_probability_matrix__(groupsrc, groupdst,
                                                      connect_dist, expand=True, hide=False)

    def __connect_by_binary_matrix__(self, groupsrc, groupdst, connect_inst,
                                     expand=True, hide=False):
        if not self.is_connect_possible(groupsrc, groupdst):
            return []

        pairs_all = np.array(np.meshgrid(groupsrc.paddr, groupdst.paddr)).transpose()
        connect_inst = np.array(connect_inst, 'bool')
        pairs_selected = pairs_all[connect_inst]

        c1 = pairs_selected[:, 0]
        c2 = pairs_selected[:, 1]

        if expand:
            self.mapping.extend(np.column_stack([c1, c2]).tolist())
        else:
            self.mapping = np.column_stack([c1, c2]).tolist()

        if not hide:
            self.add_edge(groupsrc, groupdst, arrowhead='crow', dir='both')
        
        return self.mapping

    def __connect_by_probability_matrix__(self, groupsrc, groupdst, M, expand=True, hide=False, return_inst=False):

        connect_inst = self.__instance_from_matrix_random(M)
        mapping = self.__connect_by_binary_matrix__(
            groupsrc, groupdst, connect_inst, expand=True, hide=False)

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

    def __connect_by_arbitrary_matrix__(self, groupsrc, groupdst, M, resize_method='resample', expand=True, hide=False):

        M = getattr(self, '_Mapping__resize_matrix_' + resize_method)(
            groupsrc, groupdst, M)
        connect_inst = self.__instance_from_matrix_random(M)

        return self.__connect_by_binary_matrix__(groupsrc, groupdst, connect_inst, expand=True, hide=False)

    def __connect_shuffle_all2all__(self, groupsrc, groupdst, p=0.25, expand=True, hide=False):
        """
        Connects in an all to all fashion by shuffling p*N/N 1's and 1-p*N/N 0's for every source, from the first of the source to the last of the source.
        Has the advantage that the number of connections is kept fixed, which can reduce inhomogeneity in the outputs.
        p can be an iterable, in which case the probability of each source to connect to the target i is equal to p[i].
        i.e. creates a matrix M such that P(M_{src,tgt}==1) = p[tgt] with src=1, ..., Nsrc and tgt=1, ... ,Ntgt
        """
        if not self.is_connect_possible(groupsrc, groupdst):
            return []

        if not hasattr(p, '__len__'):
            p = p * np.ones(len(groupdst))

        N = [p[i] * len(groupsrc) for i in range(len(groupdst))]
        connect_dist = np.zeros([len(groupsrc), len(groupdst)]).astype('bool')

        # It is important to know that p is intepreted as the probability to
        # connect to the target, not that the source connects to any target
        #i.e. M such that P(M_{src,tgt}==1) = p[tgt]

        #Create N[i] for each source i...
        for i in range(len(groupdst)):
            #Separate integer and factional parts
            frac, integ = np.modf(N[i])
            connect_dist[:integ, i] = True
            # Deal with fractional part by rolling a dice and checking not to
            # overrun
            if np.random.rand() < frac and integ + 1 < len(groupsrc):
                connect_dist[integ + 1, i] = True

        #... and shuffle them
        for i in range(connect_dist.shape[1]):
            np.random.shuffle(connect_dist[:, i])

        pairs_all = np.array(
            np.meshgrid(groupsrc.paddr, groupdst.paddr)).transpose()
        pairs_selected = pairs_all[connect_dist]

        c1 = pairs_selected[:, 0]
        c2 = pairs_selected[:, 1]

        if expand:
            self.mapping.extend(np.column_stack([c1, c2]))
        else:
            self.mapping = np.column_stack([c1, c2])

        if not hide:
            self.add_edge(groupsrc, groupdst, arrowhead='crow', dir='both')
        
        return self.mapping

    def add_edge(self, groupsrc, groupdst, arrowhead='normal', dir='forward'):
        scluster = pydot.Cluster(str(groupsrc.channel))
        s = pydot.Node(str(groupsrc.name))
        scluster.add_node(s)

        dcluster = pydot.Cluster(str(groupdst.channel))
        d = pydot.Node(str(groupdst.name))
        dcluster.add_node(d)

        self.graph.add_subgraph(dcluster)
        self.graph.add_subgraph(scluster)

        e = pydot.Edge(s, d, arrowhead=arrowhead, dir=dir)

        self.graph.add_edge(e)

    def save_graph(self, filename=None):
        # TODO: load graph
        if not filename:
            filename = "graph"
        f = open(filename + '.dot', 'w')
        f.write(self.graph.to_string())
        f.close()

    def clear(self):
        """
        Clear mapping list.
        """
        self.mapping = []
        self.graph = pydot.Graph('AER_connectivity')

    #def gui(self,verbose=False):
        #if len(self.mapping)==0: raise Exception('no mapping to display')
        #return MappingGui.MappingGui(self,verbose=verbose)

#Commented because of new map_api bindings
#
#    def read(self, verbose = False):
#        """
#        Donwnloads the mapping from the mapper. *** API ***
#        """
#        # TODO: subprocess doesn't work??!
# p = subprocess.Popen('ssh -C root@' + self.host + ' "cd mappinglib &&
# ./getallmappings"')
#        # p.wait()
# p = os.popen('ssh -C root@' + self.host + ' "cd mappinglib &&
# ./getallmappings"')
#
#        self.clear()
#        while True:
#            l= p.readline()
#            if l=='': break
#            x= map(lambda x:int(x), l.strip().split(' ') )
#            if len(x) == 2: self.mapping.append(x)
#        self.__toint__()
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

    def __connect_by_probability_matrix__(self, groupsrc, groupdst, M, expand=True, hide=False, return_inst=False):
        '''
        M : shape -> (len(groupsrc), len(groupdst))
        '''

        if not self.is_connect_possible(groupsrc, groupdst):
            return []

        #Create pairs of source - destinations and reshape as appropiately
        pairs = np.array(np.meshgrid(groupsrc.paddr, groupdst.paddr)
            ).transpose().reshape(-1, 2)

        p = np.column_stack([pairs, M.flatten()])

        # Transform p. Somehow this is not the best place to do this since it
        # assumes already the discretization due to the mapper firmware
        p[:, 2] = np.array(p[:, 2] * self.max_p, 'uint32')

        #Remove all zero probabilities
        zp = (p[:, 2] != 0)
        p_nonzero = p[zp, :].tolist()

        if expand:
            self.mapping.extend(p_nonzero)
        else:
            self.mapping = p_nonzero

        if not hide:
            self.add_edge(groupsrc, groupdst, arrowhead='crow', dir='both')
        
        return self.mapping

    def complete(self):
        '''
        For all connections without a probability, assume the probability is one
        '''
        m = self.mapping
        for c in m:
            if len(c) == 2:
                c.append(self.max_p)

    def prepare(self):
        super(PMapping,self).prepare()
        self.complete()
