#-----------------------------------------------------------------------------
# Purpose:
#
# Author: <authors name>
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
class ResourceManagerBase(object):
    '''
    The ResourceManagerBase class is a base class for opening, closing a resource. It is used as parent classes for the configuration, communication and mapping APIs.
    The init function takes no arguments by default
    '''
    def __init__(self):
        self._isopen = False

    def open(self):
        '''
        Opens resource
        '''
        self._isopen = True

    def close(self):
        '''
        Closes resource
        '''
        self._isopen = False

    @property
    def isopen(self):
        return self._isopen

    def __del__(self):
        if self.isopen:
            self.close()


class BatchCommunicatorBase(ResourceManagerBase):
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

    def run(self, stimulus=None, duration=None, context_manager=None):
        #IMPLEMENT
        '''
        Stimulate the neuromorphic hardware with event stream *stimulus* and monitor it for *duration* ms
        Input:
        *stimulus*: a numpy array in (addr, time) format. Type should be uint32 (shape = (-1,2)).
        *duration*: monitor duration.
        *context_manager*: context manager used to wrap the stimulate function. Useful for syncing with external process.
        Output: a numpy array in addr, time format. Type is uint32 (shape = (-1,2)).

        Usage:
        >>> events = np.array([[0, 1],[100,200]], dtype='uint32') #addr, time format
        >>> run(events, 1500)
        '''
        pass


class ContinuousCommunicatorBase(BatchCommunicatorBase):

    def __init__(self, *args, **kwargs):
        '''
        The Continuous Communicator API defined by the NCS tools.
        In addition to the Batch Communicator, this class has a mon() method which can be used to monitor data in real-time.

        Example Usage:
        >>> c = Communicator(host = 'localhost', devnum = 0)
        >>> c.run()
        >>> #Or
        >>> c = Communicator(host = 'localhost', devnum = 0)
        >>> c.open()
        >>> c.mon()
        >>> c.close()
        '''

#        Inputs:
# *host:* the hostname of the computer where the hardware is attached (default:
# localhost).
        ResourceManagerBase.__init__(self)

    def mon(self, duration):
        #IMPLEMENT
        '''
        Monitor the neuromorphic hardware for *duration* ms
        '''
        pass

    def stim(self, stimulus, duration=None, context_manager=None, **stim_kwargs):
        #IMPLEMENT
        '''
        Stimulate the neuromorphic hardware with event stream *stimulus* and monitor it for *duration* ms
        Input:
        *stimulus*: a numpy array in addr, time format. Type should be uint32 (shape = (-1,2)).
        *duration*: monitor duration.
        *context_manager*: context manager used to wrap the stimulate function. Useful for syncing with external process.
        ***stim_kwargs*: keyward arguments passed to the underlying stimulation modules

        Output: a numpy array in addr, time format. Type is uint32 (shape = (-1,2)).
        Usage:
        >>> events = np.array([0, 1],[100,200], dtype='uint32') #addr, time format
        >>> stim(events, 1500)
        '''
        pass

Communicator = BatchCommunicatorBase
