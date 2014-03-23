#-----------------------------------------------------------------------------
# Purpose:
#
# Author: <authors name>
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------

#Path for recording experimental data
from ComAPI import *
import collections 

class BatchCommunicator(BatchCommunicatorBase):
    def __init__(self):
        '''
        The BatchCommunicator API defined by the NCS tools.

        Usage:
        >>> c = Communicator(host = 'localhost', devnum = 0)
        >>> c.run()
        '''

#        Inputs:
# *host:* the hostname of the computer where the hardware is attached (default:
# localhost).
        ResourceManagerBase.__init__(self)
        RecordableCommunicatorBase.__init__(self)
        self.queue = collections.deque()

    def run(self, stimulus=None, duration=None, context_manager=None):
        '''
        Loopback API simply puts an event packet in the API, and reads it back
        '''
        self.queue.appendleft(stimulus)
        return self.queue.pop()


class ContinuousCommunicator(BatchCommunicator):

    def mon(self, duration=None):
        #IMPLEMENT
        '''
        Returns one event packet from the queue independent of the durection
        '''
        return self.queue.pop()
        

    def stim(self, stimulus, duration=None, context_manager=None, **stim_kwargs):
        #IMPLEMENT
        '''
        Adds an event packet to the queue, independent of duration and reads it.
        '''
        self.queue.appendleft(stimulus)
        return self.mon()

    def run(self, stimulus=None, duration=None, context_manager=None):
        '''
        Loopback API simply puts an event packet in the API, and reads it back
        '''
        return self.stim(stimulus)



        



Communicator = ContinuousCommunicator
