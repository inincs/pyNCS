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
import pylab

from .mapping import Mapping, PMapping


class Connection(object):
    """
    A class representing the connections between populations.
    """
    def __init__(self, popsrc, popdst, synapse,
                 fashion='one2one',
                 fashion_kwargs={},
                 connection_kwargs={},
                 append=True,
                 setup=None):
        """
        - popsrc: source population
        - popdst: destination population
        - synapse: on which synapse to connect
        - fashion: type of connectivity
        - fashion_kwargs: arguments for fashion-type connectivity
        - append: whether to append this connection to the setup mapping table
        - setup: specify setup if different from popsrc.setup
        """
        self.mapping = self._create_mapping(popsrc, popdst, synapse)
        self.mapping.connect(
                popsrc.soma,
                popdst.synapses[synapse],
                connection_kwargs = connection_kwargs,
                fashion=fashion,
                fashion_kwargs=fashion_kwargs)
        self.mapping.prepare()
        if setup is None:
            setup = popsrc.setup
        if append:
            setup.mapping.merge(self.mapping)

        self.synapse = synapse
        self.popsrc = popsrc
        self.popdst = popdst
        self.connection_kwargs = connection_kwargs

    def _create_mapping(self, popsrc, popdst, synapse):
        return Mapping(popsrc.name + '_to_' + popdst.name,
                                     'Connects %s to %s on synapse \'%s\'.' %
                                     (popsrc.name, popdst.name, synapse))

    def __len__(self):
        return len(self.mapping.mapping)

    def __repr__(self):
        return "Connection object: {0} -> {1} via {2}".format(self.popsrc.name,
                                                              self.popdst.name,
                                                              self.synapse)


    def plot(self):
        '''
        plots the connectivity
        '''
        srcindx = {j:i for i,j in enumerate(self.popsrc.soma.paddr)}
        dstindx = {j:i for i,j in enumerate(self.popdst.synapses[self.synapse].paddr)}
        conn_matrix = np.zeros((len(self.popdst.synapses[self.synapse]),
                                len(self.popsrc.soma)))
        for src, dst in np.array(self.mapping.mapping)[:,:2]:
            conn_matrix[dstindx[dst], srcindx[src]] += 1
        pylab.colorbar(pylab.pcolor(conn_matrix))
        pylab.xlabel(self.popsrc.name)
        pylab.ylabel(self.popdst.name)
        
#    def __connect_one2one__(self, popsrc, popdst, synapse, syn_ids=[0]):
#        """
#        Connects in a one to one fashion. Every source neuron connects to one
#        synapse on the destination.
#            - syn_ids: array of synapses into which to connect
#        """
#        nsrc = len(syn_ids)
#        ndst = len(popdst[0].synapses[synapse])
#        [self.mapping.connect(popsrc.soma[i::nsrc],
#                              popdst.synapses[synapse][i::ndst])
#                              for i in syn_ids]

#    def activate(self, setup=None):
#        """
#        Shortcut function to apply connections over-writing the existing ones.
#        See also NeuroSetup.Mapping.import_from_connections().
#            - setup: if None, self.setup will be used
#        """
#        if setup is None:
#            self.mapping.write()
#        else:
#            setup.mapping.import_from_connections(self)
#            setup.mapping.write()


class PConnection(Connection):
    def _create_mapping(self, popsrc, popdst, synapse):
        return PMapping(popsrc.name + '_to_' + popdst.name,
                        'Connects %s to %s on synapse \'%s\'.' %
                        (popsrc.name, popdst.name, synapse))
