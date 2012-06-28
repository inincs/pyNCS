#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Fabio Stefanini
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
import numpy as np

import pyAex
from mapping import Mapping, PMapping


class Connection(object):
    """
    A class representing the connections between populations.
    """
    def __init__(self, popsrc, popdst, synapse,
                 fashion='one2one',
                 fashion_kwargs={},
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

#        if hasattr(self, '__connect_' + fashion + '__'):
#            getattr(self, '__connect_' + fashion + '__')(popsrc, popdst, synapse,
#                                                     **fashion_kwargs)
#        else:
        self.mapping.connect(popsrc.soma,
                                 popdst.synapses[synapse],
                                 fashion=fashion,
                                 fashion_kwargs=fashion_kwargs)
        if setup is None:
            setup = popsrc.setup
        if append:
            setup.mapping.merge(self.mapping)

        self.namesrc = popsrc.name
        self.namedst = popdst.name
        self.synapse = synapse

    def _create_mapping(self, popsrc, popdst, synapse):
        return Mapping(popsrc.name + '_to_' + popdst.name,
                                     'Connects %s to %s on synapse \'%s\'.' %
                                     (popsrc.name, popdst.name, synapse))

    def __len__(self):
        return len(self.mapping.mapping)

    def __repr__(self):
        return "Connection object: {0} -> {1} via {2}".format(self.namesrc, self.namedst, self.synapse)

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
