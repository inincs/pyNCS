#-----------------------------------------------------------------------------
# Purpose:
#
# Author: <authors name>
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------

#Path for recording experimental data
import numpy, getpass, warnings
EMPTY_RAW = numpy.zeros([0,2], dtype='uint32')
USER = getpass.getuser()

class ResourceManagerBase(object):
    '''
    The ResourceManagerBase class is a base class for opening, closing a resource. It is used as parent classes for the configuration, communication and mapping APIs.
    The init function takes no arguments by default
    '''
    def __init__(self):
        self._neurosetup = None
        self._neurosetup_registered = False
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

    @property
    def neurosetup(self):
        if not self._neurosetup_registered:
            warnings.warn('NeuroSetup has not been registered')
            return None
        else:
            return self._neurosetup

    def register_neurosetup(self, neurosetup):
        '''
        Provides a link to the Neurosetup. This is useful for complex parameter
        configuration protocols requiring the sequencing and monitoring of
        address-events
        '''
        self._neurosetup_registered = True
        self._neurosetup = neurosetup

class RecordableCommunicatorBase(object):
    REC_PATH = '/tmp/exp_rec_' + USER
    REC_FN_SEQ = 'seq'
    REC_FN_MON = 'mon'
    REC_HEADER_SEQ = '# File format raw address, ISI (us)'
    REC_HEADER_MON = '# File format raw address, timestamp (us)'

    def __init__(self):
        self._rec_fns = []
        self._run_id = 0
        self.__reset()
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
        self.__reset()
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
        self.__save_rec_file(stimulus, stim_fn)
        mon_evs = self.run(stimulus = stimulus, *args, **kwargs)
        self.__save_rec_file(mon_evs, mon_fn)

        return mon_evs

    def __save_rec_file(self, ev_array, filename):
        '''
        Save data using np.save, and adds filename in self._record_fns
        '''
        self._rec_fns.append(filename+'.npy')
        numpy.save(filename, ev_array)

        self._run_id += 1

    def __gen_rec_fns(self):
        '''
        Generate filenames for recording
        '''
        #Get username


        import time
        N = self._run_id
        current_time_str = time.strftime("__" + "%d-%m-%Y-%H:%M:%S", time.localtime())
        filename = self.REC_PATH + '_{2}_{0}__run{1}'
        stim_fn = filename.format(current_time_str, N, self.REC_FN_SEQ)
        mon_fn = filename.format(current_time_str, N, self.REC_FN_MON )
        return stim_fn, mon_fn

    def __reset(self):
        self._run_id = 0
        self._rec_fns = []
        return None

    def reset(self):
        #CONVIENIENCE FUNCTION, IMPLEMENTATION NOT REQURIED
        pass
    
    def del_all(self):
        import glob, os
        fn_to_del = glob.glob(self.REC_PATH+'*')
        for f in fn_to_del:
            os.remove(f)

    def __del__(self):
        self.del_all()

    def run(self, *args, **kwargs):
        return EMPTY_RAW



