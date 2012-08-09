#-----------------------------------------------------------------------------
# Purpose:
#
# Author: <authors name>
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------

#Path for recording experimental data

import numpy
EMPTY_RAW = numpy.zeros([0,2], dtype='uint32')

class RecordableCommunicatorBase(object):
    REC_PATH = '/tmp/exp_rec'
    REC_FN_SEQ = 'seq'
    REC_FN_MON = 'mon'
    REC_HEADER_SEQ = '# File format raw address, ISI (us)'
    REC_HEADER_MON = '# File format raw address, timestamp (us)'

    def __init__(self):
        self._rec_fns = []
        self._run_id = 0
        self.reset()
        self.del_all()

    def get_exp_rec(self):
        '''
        Returns the filename where the raw stimated and monitored events are saved.
        The filenames for all experiments are returned, and run_id is reset

        Output: list of filenames containing the experiment records (stimulated and monitored events). Additionally, the experiment number is reset.
        
        '''        
        #CONVIENIENCE FUNCTION, IMPLEMENTATION NOT REQURIED
        import copy
        rec_fns = copy.copy(self._rec_fns)
        self.reset()
        return rec_fns

    def run_rec(self, stimulus = None, *args, **kwargs):
        #CONVIENIENCE FUNCTION, IMPLEMENTATION NOT REQURIED
        '''
        Stimulate the neuromorphic hardware with event stream *stimulus* and monitor it for *duration* ms.
        In addition to the run() function, this function records the experimental data in /tmp/.
        The resulting files can be obtained with get_exp_record
        This function should be overridden to avoid redundant events file creation.

        Input: (see also run() for more information)
        *stimulus*: a numpy array in (addr, time) format. Type should be uint32 (shape = (-1,2)).
        *duration*: monitor duration.
        *context_manager*: context manager used to wrap the stimulate function. Useful for syncing with external process.
        Output: a numpy array in addr, time format. Type is uint32 (shape = (-1,2)).

        Usage:
        >>> events = np.array([[0, 1],[100,200]], dtype='uint32') #addr, time format
        >>> run_record(events, 1500)
        '''
        #CONVIENIENCE FUNCTION, IMPLEMENTATION NOT REQUIRED
        stim_fn, mon_fn = self.__gen_rec_fns()
        #Save stim before in case something goes wrong        
        self.__save_rec_file(stimulus, stim_fn, header = self.REC_HEADER_SEQ)
        mon_evs = self.run(stimulus = stimulus, *args, **kwargs)
        self.__save_rec_file(mon_evs, mon_fn, header = self.REC_HEADER_MON)

        return mon_evs

    def __save_rec_file(self, ev_array, filename, *args, **kwargs):
        '''
        Save data using np.savetxt, and adds filename in self._record_fns
        '''
        import numpy
        if int(numpy.__version__.split('.')[0])==1 and int(numpy.__version__.split('.')[1])<=6:
            kwargs.pop('header')
        self._rec_fns.append(filename)
        numpy.savetxt(filename, ev_array, delimiter = ' ', newline = '\n', fmt = ['%u', '%u'], *args, **kwargs)
        self._run_id += 1

    def __gen_rec_fns(self):
        '''
        Generate filenames for recording
        '''
        import time
        N = self._run_id
        current_time_str = time.strftime("__" + "%d-%m-%Y-%H:%M:%S", time.localtime())
        filename = self.REC_PATH + '_{2}_{0}__run{1}.dat'
        stim_fn = filename.format(current_time_str, N, self.REC_FN_SEQ)
        mon_fn = filename.format(current_time_str, N, self.REC_FN_MON )
        return stim_fn, mon_fn

    def reset(self):
        self._run_id = 0
        self._rec_fns = []
        return None
    
    def del_all(self):
        import glob, os
        fn_to_del = glob.glob(self.REC_PATH+'*')
        for f in fn_to_del:
            os.remove(f)

    def __del__(self):
        self.del_all()

    def run(self, *args, **kwargs):
        return EMPTY_RAW


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


class BatchCommunicatorBase(ResourceManagerBase,RecordableCommunicatorBase):
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
        return EMPTY_RAW




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
        BatchCommunicatorBase.__init__(self)

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
        return EMPTY_RAW

Communicator = BatchCommunicatorBase
