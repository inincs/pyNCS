import numpy as np

class Events(object):
    
    NF = 2 #Number of dimensions (address, timestamp)

    def __init__(self, ev=None, atype='p', isISI=False):
        '''
        Class to hold event data related data types.
        '''
        self.isISI = isISI
        assert isinstance(atype, str)
        __atype = atype[0].lower()

        if __atype == 'p':
            __dtype = np.dtype([('tm', 'uint32'),
                              ('ad', 'uint32')])

        if __atype == 'l':
            __dtype = np.dtype([('tm', 'float'),
                              ('ad', 'float')])

        self.__atype = __atype
        self.__dtype = __dtype

        if isinstance(ev, events):
            self.__data = ev.__data.copy()
        elif ev != None:
            ev = np.array(ev)
            if ev.shape[1] == self.NF:
                self.set_data(ev[:, 0], ev[:, 1])
            elif ev.shape[0] == self.NF:
                self.set_data(ev[0, :], ev[1, :])
            else:
                raise TypeError
        else:
            self.__data = np.zeros([0], self.dtype)

    @property
    def data(self):
        return self.__data

    @property
    def atype(self):
        return self.__atype

    @property
    def dtype(self):
        return self.__dtype

    def __len__(self):
        return self.get_nev()

    def __add__(self, other):
        self.add_adtm(other.ad, other.tm)

    def __repr__(self):
        return self.get_adtm().__repr__()

    def get_nev(self):
        return len(self.get_tm())

    @property
    def nev(self):
        return self.get_nev()

    def get_tdur(self):
        if self.nev == 0:
            return 0
        if self.isISI:
            return self.tm.sum()
        else:
            if self.get_nev() > 0:
                return self.tm[-1] - self.tm[0]

    @property
    def tdur(self):
        return self.get_tdur()

    def get_t_stop(self):
        if self.isISI:
            return self.tm.sum()
        else:
            if self.get_nev() > 0:
                return self.tm[-1]

    @property
    def t_stop(self):
        return self.get_t_stop()

    def get_ad(self):
        return self.__data['ad']

    def set_ad(self, ad):
        self.__data['ad'] = ad
        return self.get_ad()

    @property
    def ad(self):
        return self.get_ad()

    #@ad.setter
    #def ad(self, value):
    #    return self.set_ad(value)

    def get_tm(self):
        return self.__data['tm']

    def set_tm(self, tm):
        self.__data['tm'] = tm
        return self.get_tm()

    @property
    def tm(self):
        return self.get_tm()

    #@tm.setter
    #def tm(self, value):
    #    return self.set_tm(value)

    def set_data(self, ad, tm):
        assert len(ad) == len(tm), "addresses and timestamps lengths are incompatible %d %d" % (len(ad), len(tm))
        self.__data = np.zeros(len(ad), self.dtype)
        self.__data['ad'] = ad.astype(self.dtype['ad'])
        self.__data['tm'] = tm.astype(self.dtype['tm'])

    def add_adtmev(self, ev):
        if not isinstance(ev, np.ndarray):
            ev = np.array(ev)

        if len(ev.shape) != self.NF:
            ev = ev.reshape(-1, 2).astype(self.dtype)

        ad = np.concatenate([self.ad, ev[:, 0]])
        tm = np.concatenate([self.tm, ev[:, 1]])
        self.set_data(ad, tm)

    def add_adtm(self, ad, tm):
        if not isinstance(ad, np.ndarray):
            ad = np.array(ad)

        if len(ad.shape) != 1:
            ad = ad.reshape(-1, 1)

        if not isinstance(tm, np.ndarray):
            tm = np.array(tm)

        if len(tm.shape) != 1:
            tm = tm.reshape(-1, 1)

        assert tm.shape == ad.shape

        ad = np.concatenate([self.ad, ad])
        tm = np.concatenate([self.tm, tm])
        self.set_data(ad, tm)

    def get_tmad(self):
        return np.array([self.tm, self.ad])

    def get_adtm(self):
        return np.array([self.ad, self.tm])

    def get_tmadev(self):
        return self.get_tmad().transpose()

    def get_adtmev(self):
        return self.get_adtm().transpose()

    def normalize_tm(self, t0=0.):
        if not self.isISI:
            t_start = self.get_tm().min()
            self.set_tm(self.get_tm() - t_start + t0)
        else:
            t = self.get_tm()
            t[0] = t0
            self.set_tm(t)

    def get_adisi(self):

        if self.isISI:
            return self.get_adtm()
        else:
            if self.nev > 0:
                tm = np.concatenate([np.array([self.tm[0]]), np.diff(self.tm)])
                return np.array([self.ad, tm])
            else:
                return np.zeros([2,0])

    def get_isiad(self):

        if self.isISI:
            return self.get_tmad()
        else:
            if self.nev > 0:
                tm = np.concatenate([np.array([self.tm[0]]), np.diff(self.tm)])
                return np.array([tm, self.ad])
            else:
                return np.zeros([2,0])

    def set_isi(self):
        if self.isISI:
            pass
        else:
            evs = self.get_adisi()
            self.set_data(evs[0],evs[1])
            self.isISI = True
            

    def set_abs_tm(self):
        '''
        Transform ISI timestamps into absolute time 
        '''
        if self.isISI:
            self.ad, self.set_tm(np.cumsum(self.get_tm()))
            self.isISI = False
        else:
            pass

    def filter_by_mapping(self, mapping):
        """
        Map the events, given a mapping dictionary like:
        map[src]=[target1,target2,...,targetn],
        """
        #The following is optimized for performance
        wasISI = False
        if self.isISI:
            wasISI = True
            self.set_abs_tm()
        evs = self.get_adtmev()
        ff = lambda x: x[0] in mapping
        filt = filter(ff, evs) #keep only addresses that are mapped
        if len(filt) > 0:
            evs_filt = events(filt, self.atype).get_adtmev()
            #Get mapped addresses
            #list(chain(* concatenates lists of lists for low cost (i.e. O(n))
            m_ad = np.array(list(itertools.chain(*map(mapping.get, evs_filt[:, 0]))), self.dtype['ad'])
            #Get the respective timestamps
            m_tm = np.array(list(itertools.chain(*map(lambda x: len(mapping.get(x[0])) * [x[1]], evs_filt))), self.dtype['tm'])

            self.set_data(m_ad,m_tm)
        else:
            self.empty()
        if wasISI:
            self.set_isi()

    def empty(self):
        self.__data = np.zeros([0], self.dtype)

    def iter_by_timeslice(self, tm):
        import bisect
        #ISI -> Cumulative
        if self.isISI:
            sum_tm = np.cumsum(self.get_tm())
        else:
            sum_tm = self.get_tm()
        #No in-place change
        id_start = 0
        t = 0

        #better to recycle an object rather than creating new (slow)
        evs = events(isISI=self.isISI)

        while id_start < len(sum_tm):
            t += tm
            id_stop = bisect.bisect_right(sum_tm, t, lo=id_start)
            evs.__data = self.__data[id_start:id_stop]
            id_start = id_stop
            rest = tm - evs.get_tdur()
            print(tm, evs.get_tdur())

            if not evs.get_nev() > 0:
                rest = 0

            t -= rest
            #ISI or not is determined
            yield evs, tm - rest

    def __iter__(self):
        for ev in self.get_tmadev():
            yield ev[0],ev[1]

    def sort(self):
        '''
        Sort events by addresses and timestamps (in this order).
        '''
        self.set_abs_tm()
        self.__data.sort(order=['ad', 'tm'])

    def demultiplex(self):
        '''
        Generates a dictionary with addesses as keys and a list of timestamps as values.
        Used internally for generating SpikeLists
        '''
        evs = events(ev=self)
        evs.sort()
        ads = np.unique(evs.ad)

        d = dict()
        k_start = 0
        k_stop = np.searchsorted(self.ad, ads, side='right')

        for i, a in enumerate(ads):
            d[a] = self.tm[k_start:k_stop[i]]
            k_start = k_stop[i]

        return d


