#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Emre Neftci
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
from __future__ import absolute_import
from .STas import events, channelEvents, RawOutput, channelAddressing, addrSpec
from .STas import setDefaultMonChannelAddress, setDefaultSeqChannelAddress,\
                 getDefaultMonChannelAddress, getDefaultSeqChannelAddress
from .STsl import STCreate
from .spikes import SpikeList, SpikeTrain, merge, merge_spikelists, \
                   merge_sequencers
from . import pyST_globals
import numpy as np


def loadtxt(filename, comments, format='t'):
    data = np.loadtxt(filename, comments=comments)
    if format == 't':
        data = np.fliplr(data)
    ev = events(data, atype='l')
    id_list = np.unique(ev.get_ad())
    return SpikeList(ev.get_adtmev(), id_list)
