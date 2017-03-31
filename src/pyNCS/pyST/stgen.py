#-----------------------------------------------------------------------------
# Purpose:
#
# Adapted from NeuroTools
#
# Licence : GPLv2
#-----------------------------------------------------------------------------
"""
NeuroTools.stgen
================

A collection of tools for stochastic process generation.


Classes
-------

StGen - Object to generate stochastic processes of various kinds
        and return them as SpikeTrain or AnalogSignal objects.


Functions
---------

shotnoise_fromspikes - Convolves the provided spike train with shot decaying exponential.

gamma_hazard - Compute the hazard function for a gamma process with parameters a,b.
"""
from __future__ import absolute_import

from .spikes import SpikeTrain
from numpy import array, log
import numpy


def gamma_hazard(x, a, b, dt=1e-4):
    """
    Compute the hazard function for a gamma process with parameters a,b
    where a and b are the parameters of the gamma PDF:
    y(t) = x^(a-1) \exp(-x/b) / (\Gamma(a)*b^a)

    Inputs:
        x   - in units of seconds
        a   - dimensionless
        b   - in units of seconds

    See also:
        inh_gamma_generator
    """

    # This algorithm is presently not used by
    # inh_gamma_generator as it has numerical problems
    # Try:
    # plot(stgen.gamma_hazard(arange(0,1000.0,0.1),10.0,1.0/50.0))
    # and look for the kinks.

    from scipy.special import gammaincc
    Hpre = -log(gammaincc(a, (x - dt) / b))
    Hpost = -log(gammaincc(a, (x + dt) / b))
    val = 0.5 * (Hpost - Hpre) / dt

    if isinstance(val, numpy.ndarray):
        val[numpy.isnan(val)] = 1.0 / b
        return val
    elif numpy.isnan(val):
        return 1.0 / b
    else:
        return val


#def gamma_hazard_rpy(x, a, b, dt=1e-4):
#    """
#    Compute the hazard function for a gamma process with parameters a,b
#    where a and b are the parameters of the gamma PDF:
#    y(t) = x^(a-1) \exp(-x/b) / (\Gamma(a)*b^a)
#
#    Inputs:
#        x   - in units of seconds
#        a   - dimensionless
#        b   - in units of seconds
#
#    See also:
#        inh_gamma_generator
#
#    """
#
#    # Used by inh_gamma_generator
#
#    # Ideally, I would like to see an implementation which does not depend on
#    # RPy
#    # but the gamma_hazard_scipy above using scipy exhibits numerical problems,
#    # as it does not
#    # support directly returning the log.
#
#    from rpy import r
#
#    # scipy.special.gammaincc has numerical problems
#    #Hpre = -log(scipy.special.gammaincc(a,(x-dt)/b))
#    #Hpost = -log(scipy.special.gammaincc(a,(x+dt)/b))
#
#    # reverting to the good old r.pgamma
#    Hpre = -r.pgamma(x - dt, shape=a, scale=b, lower=False, log=True)
#    Hpost = -r.pgamma(x + dt, shape=a, scale=b, lower=False, log=True)
#    val = 0.5 * (Hpost - Hpre) / dt
#
#    return val


class StGen:

    def __init__(self, rng=None, seed=None):
        """
        Stochastic Process Generator
        ============================

        Object to generate stochastic processes of various kinds
        and return them as SpikeTrain or AnalogSignal objects.


        Inputs:
            rng - The random number generator state object (optional). Can be None, or
                  a numpy.random.RandomState object, or an object with the same
                  interface.

            seed - A seed for the rng (optional).

        If rng is not None, the provided rng will be used to generate random numbers,
        otherwise StGen will create its own random number generator.
        If a seed is provided, it is passed to rng.seed(seed)

        Examples:
            >> x = StGen()



        StGen Methods:

        Spiking point processes:
        ------------------------

        poisson_generator - homogeneous Poisson process
        inh_poisson_generator - inhomogeneous Poisson process (time varying rate)
        inh_gamma_generator - inhomogeneous Gamma process (time varying a,b)
        inh_adaptingmarkov_generator - inhomogeneous adapting markov process (time varying)
        inh_2Dadaptingmarkov_generator - inhomogeneous adapting and
                                         refractory markov process (time varying)

        Continuous time processes:
        --------------------------

        OU_generator - Ohrnstein-Uhlenbeck process


        See also:
          shotnoise_fromspikes

        """

        if rng == None:
            self.rng = numpy.random.RandomState()
        else:
            self.rng = rng

        if seed != None:
            self.rng.seed(seed)
        self.rpy_checked = False

    def seed(self, seed):
        """ seed the gsl rng with a given seed """
        self.rng.seed(seed)

    def regular_gaussian_generator(self, rate, phase=0.0, scale=5., t_start=0.0, t_stop=1000.0, array=False):
        """
        Returns a SpikeTrain whose spikes are regularly spaced, but jittered
        according to a Gaussian distribution around the spiking period, with
        the given rate (Hz) and stopping time t_stop (milliseconds).

        Note: t_start is always 0.0, thus all realizations are as if they
        spiked at t=0.0, though this spike is not included in the SpikeList.

        Inputs:
            rate    - the rate of the discharge (in Hz)
            t_start - the beginning of the SpikeTrain (in ms)
            t_stop  - the end of the SpikeTrain (in ms)
            phase   - Offset the spiketrain by this number (in ms.)
            scale   - width of the Gaussian distribution placed at the regular
                      spike times, according to which the spike will be drawn (in ms)
            array   - if True, a numpy array of sorted spikes is returned,
                      rather than a SpikeTrain object.

        Examples:
            >> regular_generator(50, 0, 1000)
            >> regular_generator(20, 5000, 10000, array=True)

        See also:
            inh_poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator
        """

        import scipy.stats

        n = (t_stop - t_start) / 1000.0 * rate
        t_bin = 1
        time_axis = numpy.arange(t_start, t_stop, t_bin)

        if n > 0:
            t_peaks = numpy.arange(
                t_start, t_stop, 1000. / rate, dtype='float') + phase
        else:
            t_peaks = numpy.array([])

        #Some spikes will fall off the starting edge. Get rid of them
        t_peaks = t_peaks[t_peaks > 0.]

        spike_times = numpy.array(
            [scipy.stats.norm.rvs(t, scale) for t in t_peaks])

        if not array:
            return SpikeTrain(spike_times, t_start=t_start, t_stop=t_stop)
        else:
            return numpy.sort(spike_times)

    def regular_generator(self, rate, phase=0.0, jitter=True, t_start=0.0, t_stop=1000.0, array=False):
        """
        Returns a SpikeTrain whose spikes are regularly spaced
        with the given rate (Hz) and stopping time t_stop (milliseconds).

        Note: t_start is always 0.0, thus all realizations are as if
        they spiked at t=0.0, though this spike is not included in the SpikeList.

        Inputs:
            rate    - the rate of the discharge (in Hz)
            t_start - the beginning of the SpikeTrain (in ms)
            phase   - Offset the spiketrain by this number (in ms.)
            jitter  - whether the spiketrain should be jittered by an amount numpy.random.rand()/rate
            t_stop  - the end of the SpikeTrain (in ms)
            array   - if True, a numpy array of sorted spikes is returned,
                      rather than a SpikeTrain object.

        Examples:
            >> regular_generator(50, 0, 1000)
            >> regular_generator(20, 5000, 10000, array=True)

        See also:
            inh_poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator
        """

        # less wasteful than double length method above
        n = (t_stop - t_start) / 1000.0 * rate

        #jitter

        if n > 0:
            spikes = numpy.arange(
                t_start, t_stop, 1000. / rate, dtype='float') + phase
        else:
            spikes = numpy.array([])

        if jitter:
            spikes += numpy.random.rand() * 1000. / rate
            #Remove any spikes that extend beyond t_stop
            spikes = spikes[spikes<t_stop]

        if not array:
            spikes = SpikeTrain(spikes, t_start=t_start, t_stop=t_stop)

        return spikes

    def poisson_generator(self, rate, t_start=0.0, t_stop=1000.0, array=False, debug=False, refractory=0):
        """
        Returns a SpikeTrain whose spikes are a realization of a Poisson process
        with the given rate (Hz) and stopping time t_stop (milliseconds).

        Note: t_start is always 0.0, thus all realizations are as if
        they spiked at t=0.0, though this spike is not included in the SpikeList.

        Inputs:
            rate    - the rate of the discharge (in Hz)
            t_start - the beginning of the SpikeTrain (in ms)
            t_stop  - the end of the SpikeTrain (in ms)
            array   - if True, a numpy array of sorted spikes is returned,
                      rather than a SpikeTrain object.

        Examples:
            >> gen.poisson_generator(50, 0, 1000)
            >> gen.poisson_generator(20, 5000, 10000, array=True)

        See also:
            inh_poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator
        """

        #number = int((t_stop-t_start)/1000.0*2.0*rate)

        # less wasteful than double length method above
        n = (t_stop - t_start) / 1000.0 * rate
        number = int(numpy.ceil(n + 3 * numpy.sqrt(n)))
        if number < 100:
            number = min(5 + int(numpy.ceil(2 * n)), 100)

        if number > 0:
            isi = self.rng.exponential(1.0 / rate, number) * 1000.0

            if number > 1:
                spikes = numpy.add.accumulate(isi)
            else:
                spikes = isi
        else:
            spikes = numpy.array([])

        spikes += t_start
        i = numpy.searchsorted(spikes, t_stop)

        extra_spikes = []
        if i == len(spikes):
            # ISI buf overrun

            t_last = spikes[-1] + self.rng.exponential(1.0 /
                 rate, 1)[0] * 1000.0

            while (t_last < t_stop):
                extra_spikes.append(t_last)
                t_last += self.rng.exponential(1.0 / rate, 1)[0] * 1000.0

            spikes = numpy.concatenate((spikes, extra_spikes))

            if debug:
                print("ISI buf overrun handled." +
                      "len(spikes)={0}, len(extra_spikes)={1}".format(
                          len(spikes), len(extra_spikes)))

        else:
            spikes = numpy.resize(spikes, (i,))

        if not array:
            spikes = SpikeTrain(spikes, t_start=t_start, t_stop=t_stop)

        if debug:
            return spikes, extra_spikes
        else:
            return spikes

    def inh_poisson_generator(self, rate, t, t_stop, base_generator=None, array=False, **base_generator_kwargs):
        """
        Returns a SpikeList whose spikes are a realization of an inhomogeneous
        poisson process (dynamic rate). The implementation uses the thinning
        method, as presented in the references.

        Inputs:
            rate   - an array of the rates (Hz) where rate[i] is active on interval
                     [t[i],t[i+1]]
            t      - an array specifying the time bins (in milliseconds) at which to
                     specify the rate
            t_stop - length of time to simulate process (in ms)
            array  - if True, a numpy array of sorted spikes is returned,
                     rather than a SpikeList object.

        Note:
            t_start=t[0]

        References:

        Eilif Muller, Lars Buesing, Johannes Schemmel, and Karlheinz Meier
        Spike-Frequency Adapting Neural Ensembles: Beyond Mean Adaptation and Renewal Theories
        Neural Comput. 2007 19: 2958-3010.

        Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.

        Examples:
            >> time = arange(0,1000)
            >> stgen.inh_poisson_generator(time,sin(time), 1000)

        See also:
            poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator
        """

        if base_generator == None:
            base_generator = self.poisson_generator

        if numpy.shape(t) != numpy.shape(rate):
            raise ValueError(
                'shape mismatch: t,rate must be of the same shape')

        # get max rate and generate poisson process to be thinned
        rmax = numpy.max(rate)
        ps = base_generator(rmax, t_start=t[0], t_stop=t_stop,
             array=True, **base_generator_kwargs)

        # return empty if no spikes
        if len(ps) == 0:
            if array:
                return numpy.array([])
            else:
                return SpikeTrain(numpy.array([]), t_start=t[0], t_stop=t_stop)

        # gen uniform rand on 0,1 for each spike
        rn = numpy.array(self.rng.uniform(0, 1, len(ps)))

        # instantaneous rate for each spike

        idx = numpy.searchsorted(t, ps) - 1
        spike_rate = rate[idx]

        # thin and return spikes
        spike_train = ps[rn < spike_rate / rmax]

        if array:
            return spike_train

        return SpikeTrain(spike_train, t_start=t[0], t_stop=t_stop)

    def _inh_gamma_generator_python(self, a, b, t, t_stop, array=False):
        """
        Returns a SpikeList whose spikes are a realization of an inhomogeneous gamma process
        (dynamic rate). The implementation uses the thinning method, as presented in the
        references.

        Inputs:
            a,b    - arrays of the parameters of the gamma PDF where a[i] and b[i]
                     will be active on interval [t[i],t[i+1]]
            t      - an array specifying the time bins (in milliseconds) at which to
                     specify the rate
            t_stop - length of time to simulate process (in ms)
            array  - if True, a numpy array of sorted spikes is returned,
                     rather than a SpikeList object.

        Note:
            t_start=t[0]
            a is a dimensionless quantity > 0, but typically on the order of 2-10.
            a = 1 results in a poisson process.
            b is assumed to be in units of 1/Hz (seconds).

        References:

        Eilif Muller, Lars Buesing, Johannes Schemmel, and Karlheinz Meier
        Spike-Frequency Adapting Neural Ensembles: Beyond Mean Adaptation and Renewal Theories
        Neural Comput. 2007 19: 2958-3010.

        Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.

        Examples:
            See source:trunk/examples/stgen/inh_gamma_psth.py

        See also:
            inh_poisson_generator, gamma_hazard
        """

        from numpy import shape

        if shape(t) != shape(a) or shape(a) != shape(b):
            raise ValueError('shape mismatch: t,a,b must be of the same shape')

        # get max rate and generate poisson process to be thinned
        rmax = numpy.max(1.0 / b)
        ps = self.poisson_generator(
            rmax, t_start=t[0], t_stop=t_stop, array=True)

        # return empty if no spikes
        if len(ps) == 0:
            if array:
                return numpy.array([])
            else:
                return SpikeTrain(numpy.array([]), t_start=t[0], t_stop=t_stop)

        # gen uniform rand on 0,1 for each spike
        rn = numpy.array(self.rng.uniform(0, 1, len(ps)))

        # instantaneous a,b for each spike

        idx = numpy.searchsorted(t, ps) - 1
        spike_a = a[idx]
        spike_b = b[idx]

        keep = numpy.zeros(shape(ps), bool)

        # thin spikes

        i = 0
        t_last = 0.0
        t_i = 0

        while(i < len(ps)):
            # find index in "t" time
            t_i = numpy.searchsorted(t[t_i:], ps[i], 'right') - 1 + t_i
            if rn[i] < gamma_hazard((ps[i] - t_last) / 1000.0, a[t_i], b[t_i]) / rmax:
                # keep spike
                t_last = ps[i]
                keep[i] = True
            i += 1

        spike_train = ps[keep]

        if array:
            return spike_train

        return SpikeTrain(spike_train, t_start=t[0], t_stop=t_stop)

    # use slow python implementation for the time being
    # TODO: provide optimized C/weave implementation if possible
    def inh_gamma_generator(self, a, b, t, t_stop, array=False):
        """
        Returns a SpikeList whose spikes are a realization of an inhomogeneous gamma process
        (dynamic rate). The implementation uses the thinning method, as presented in the
        references.

        Inputs:
            a,b    - arrays of the parameters of the gamma PDF where a[i] and b[i]
                     will be active on interval [t[i],t[i+1]]
            t      - an array specifying the time bins (in milliseconds) at which to
                     specify the rate
            t_stop - length of time to simulate process (in ms)
            array  - if True, a numpy array of sorted spikes is returned,
                     rather than a SpikeList object.

        Note:
            t_start=t[0]
            a is a dimensionless quantity > 0, but typically on the order of 2-10.
            a = 1 results in a poisson process.
            b is assumed to be in units of 1/Hz (seconds).

        References:

        Eilif Muller, Lars Buesing, Johannes Schemmel, and Karlheinz Meier
        Spike-Frequency Adapting Neural Ensembles: Beyond Mean Adaptation and Renewal Theories
        Neural Comput. 2007 19: 2958-3010.

        Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.

        Examples:
            See source:trunk/examples/stgen/inh_gamma_psth.py

        See also:
            inh_poisson_generator, gamma_hazard
        """

        return self._inh_gamma_generator_python(a, b, t, t_stop, array)

    def _inh_adaptingmarkov_generator_python(self, a, bq, tau, t, t_stop, array=False):

        """
        Returns a SpikeList whose spikes are an inhomogeneous
        realization (dynamic rate) of the so-called adapting markov
        process (see references). The implementation uses the thinning
        method, as presented in the references.

        This is the 1d implementation, with no relative refractoriness.
        For the 2d implementation with relative refractoriness,
        see the inh_2dadaptingmarkov_generator.

        Inputs:
            a,bq    - arrays of the parameters of the hazard function where a[i] and bq[i]
                     will be active on interval [t[i],t[i+1]]
            tau    - the time constant of adaptation (in milliseconds).
            t      - an array specifying the time bins (in milliseconds) at which to
                     specify the rate
            t_stop - length of time to simulate process (in ms)
            array  - if True, a numpy array of sorted spikes is returned,
                     rather than a SpikeList object.

        Note:
            - t_start=t[0]

            - a is in units of Hz.  Typical values are available
              in Fig. 1 of Muller et al 2007, a~5-80Hz (low to high stimulus)

            - bq here is taken to be the quantity b*q_s in Muller et al 2007, is thus
              dimensionless, and has typical values bq~3.0-1.0 (low to high stimulus)

            - tau_s has typical values on the order of 100 ms


        References:

        Eilif Muller, Lars Buesing, Johannes Schemmel, and Karlheinz Meier
        Spike-Frequency Adapting Neural Ensembles: Beyond Mean Adaptation and Renewal Theories
        Neural Comput. 2007 19: 2958-3010.

        Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.

        Examples:
            See source:trunk/examples/stgen/inh_2Dmarkov_psth.py


        See also:
            inh_poisson_generator, inh_gamma_generator, inh_2dadaptingmarkov_generator

        """

        from numpy import shape

        if shape(t) != shape(a) or shape(a) != shape(bq):
            raise ValueError('shape mismatch: t,a,b must be of the same shape')

        # get max rate and generate poisson process to be thinned
        rmax = numpy.max(a)
        ps = self.poisson_generator(
            rmax, t_start=t[0], t_stop=t_stop, array=True)

        isi = numpy.zeros_like(ps)
        isi[1:] = ps[1:] - ps[:-1]
        isi[0] = ps[0]  # -0.0 # assume spike at 0.0

        # return empty if no spikes
        if len(ps) == 0:
            return SpikeTrain(numpy.array([]), t_start=t[0], t_stop=t_stop)

        # gen uniform rand on 0,1 for each spike
        rn = numpy.array(self.rng.uniform(0, 1, len(ps)))

        # instantaneous a,bq for each spike

        idx = numpy.searchsorted(t, ps) - 1
        spike_a = a[idx]
        spike_bq = bq[idx]

        keep = numpy.zeros(shape(ps), bool)

        # thin spikes

        i = 0
        t_last = 0.0
        t_i = 0
        # initial adaptation state is unadapted, i.e. large t_s
        t_s = 1000 * tau

        while(i < len(ps)):
            # find index in "t" time, without searching whole array each time
            t_i = numpy.searchsorted(t[t_i:], ps[i], 'right') - 1 + t_i

            # evolve adaptation state
            t_s += isi[i]

            if rn[i] < a[t_i] * numpy.exp(-bq[t_i] * numpy.exp(-t_s / tau)) / rmax:
                # keep spike
                keep[i] = True
                # remap t_s state
                t_s = -tau * numpy.log(numpy.exp(-t_s / tau) + 1)
            i += 1

        spike_train = ps[keep]

        if array:
            return spike_train

        return SpikeTrain(spike_train, t_start=t[0], t_stop=t_stop)

    # use slow python implementation for the time being
    # TODO: provide optimized C/weave implementation if possible
    inh_adaptingmarkov_generator = _inh_adaptingmarkov_generator_python

    def _inh_2Dadaptingmarkov_generator_python(self, a, bq, tau_s, tau_r, qrqs, t, t_stop, array=False):

        """
        Returns a SpikeList whose spikes are an inhomogeneous
        realization (dynamic rate) of the so-called 2D adapting markov
        process (see references).  2D implies the process has two
        states, an adaptation state, and a refractory state, both of
        which affect its probability to spike.  The implementation
        uses the thinning method, as presented in the references.

        For the 1d implementation, with no relative refractoriness,
        see the inh_adaptingmarkov_generator.

        Inputs:
            a,bq    - arrays of the parameters of the hazard function where a[i] and bq[i]
                     will be active on interval [t[i],t[i+1]]
            tau_s    - the time constant of adaptation (in milliseconds).
            tau_r    - the time constant of refractoriness (in milliseconds).
            qrqs     - the ratio of refractoriness conductance to adaptation conductance.
                       typically on the order of 200.
            t      - an array specifying the time bins (in milliseconds) at which to
                     specify the rate
            t_stop - length of time to simulate process (in ms)
            array  - if True, a numpy array of sorted spikes is returned,
                     rather than a SpikeList object.

        Note:
            - t_start=t[0]

            - a is in units of Hz.  Typical values are available
              in Fig. 1 of Muller et al 2007, a~5-80Hz (low to high stimulus)

            - bq here is taken to be the quantity b*q_s in Muller et al 2007, is thus
              dimensionless, and has typical values bq~3.0-1.0 (low to high stimulus)

            - qrqs is the quantity q_r/q_s in Muller et al 2007,
              where a value of qrqs = 3124.0nS/14.48nS = 221.96 was used.

            - tau_s has typical values on the order of 100 ms
            - tau_r has typical values on the order of 2 ms


        References:

        Eilif Muller, Lars Buesing, Johannes Schemmel, and Karlheinz Meier
        Spike-Frequency Adapting Neural Ensembles: Beyond Mean Adaptation and Renewal Theories
        Neural Comput. 2007 19: 2958-3010.

        Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.

        Examples:
            See source:trunk/examples/stgen/inh_2Dmarkov_psth.py

        See also:
            inh_poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator

        """

        from numpy import shape

        if shape(t) != shape(a) or shape(a) != shape(bq):
            raise ValueError('shape mismatch: t,a,b must be of the same shape')

        # get max rate and generate poisson process to be thinned
        rmax = numpy.max(a)
        ps = self.poisson_generator(
            rmax, t_start=t[0], t_stop=t_stop, array=True)

        isi = numpy.zeros_like(ps)
        isi[1:] = ps[1:] - ps[:-1]
        isi[0] = ps[0]  # -0.0 # assume spike at 0.0

        # return empty if no spikes
        if len(ps) == 0:
            return SpikeTrain(numpy.array([]), t_start=t[0], t_stop=t_stop)

        # gen uniform rand on 0,1 for each spike
        rn = numpy.array(self.rng.uniform(0, 1, len(ps)))

        # instantaneous a,bq for each spike

        idx = numpy.searchsorted(t, ps) - 1
        spike_a = a[idx]
        spike_bq = bq[idx]

        keep = numpy.zeros(shape(ps), bool)

        # thin spikes

        i = 0
        t_last = 0.0
        t_i = 0
        # initial adaptation state is unadapted, i.e. large t_s
        t_s = 1000 * tau_s
        t_r = 1000 * tau_s

        while(i < len(ps)):
            # find index in "t" time, without searching whole array each time
            t_i = numpy.searchsorted(t[t_i:], ps[i], 'right') - 1 + t_i

            # evolve adaptation state
            t_s += isi[i]
            t_r += isi[i]

            if rn[i] < a[t_i] * numpy.exp(-bq[t_i] * (numpy.exp(-t_s / tau_s) + qrqs * numpy.exp(-t_r / tau_r))) / rmax:
                # keep spike
                keep[i] = True
                # remap t_s state
                t_s = -tau_s * numpy.log(numpy.exp(-t_s / tau_s) + 1)
                t_r = -tau_r * numpy.log(numpy.exp(-t_r / tau_r) + 1)
            i += 1

        spike_train = ps[keep]

        if array:
            return spike_train

        return SpikeTrain(spike_train, t_start=t[0], t_stop=t_stop)

    # use slow python implementation for the time being
    # TODO: provide optimized C/weave implementation if possible
    inh_2Dadaptingmarkov_generator = _inh_2Dadaptingmarkov_generator_python

    def _OU_generator_python(self, dt, tau, sigma, y0, t_start=0.0, t_stop=1000.0, array=True, time_it=False):
        """
        Generates an Orstein Ulbeck process using the forward euler method. The function returns
        an AnalogSignal object.

        Inputs:
            dt      - the time resolution in milliseconds of th signal
            tau     - the correlation time in milliseconds
            sigma   - std dev of the process
            y0      - initial value of the process, at t_start
            t_start - start time in milliseconds
            t_stop  - end time in milliseconds
            array   - if True, the functions returns the tuple (y,t)
                      where y and t are the OU signal and the time bins, respectively,
                      and are both numpy arrays.

        Examples:
            >> stgen.OU_generator(0.1, 2, 3, 0, 0, 10000)

        See also:
            OU_generator_weave1
        """

        import time

        if time_it:
            t1 = time.time()

        t = numpy.arange(t_start, t_stop, dt)
        N = len(t)
        y = numpy.zeros(N, float)
        gauss = self.rng.standard_normal(N - 1)
        y[0] = y0
        fac = dt / tau
        noise = numpy.sqrt(2 * fac) * sigma

        # python loop... bad+slow!
        for i in xrange(1, N):
            y[i] = y[i - 1] + fac * (y0 - y[i - 1]) + noise * gauss[i - 1]

        if time_it:
            print(time.time()-1)

        if array:
            return (y, t)
        else:
            raise NotImplementedError()

    # use slow python implementation for the time being
    # TODO: provide optimized C/weave implementation if possible

    def _OU_generator_python2(self, dt, tau, sigma, y0, t_start=0.0, t_stop=1000.0, array=False, time_it=False):
        """
        Generates an Orstein Ulbeck process using the forward euler method. The function returns
        an AnalogSignal object.

        Inputs:
            dt      - the time resolution in milliseconds of th signal
            tau     - the correlation time in milliseconds
            sigma   - std dev of the process
            y0      - initial value of the process, at t_start
            t_start - start time in milliseconds
            t_stop  - end time in milliseconds
            array   - if True, the functions returns the tuple (y,t)
                      where y and t are the OU signal and the time bins, respectively,
                      and are both numpy arrays.

        Examples:
            >> stgen.OU_generator(0.1, 2, 3, 0, 0, 10000)

        See also:
            OU_generator_weave1
        """

        import time

        if time_it:
            t1 = time.time()

        t = numpy.arange(t_start, t_stop, dt)
        N = len(t)
        y = numpy.zeros(N, float)
        y[0] = y0
        fac = dt / tau
        gauss = fac * y0 + numpy.sqrt(
            2 * fac) * sigma * self.rng.standard_normal(N - 1)
        mfac = 1 - fac

        # python loop... bad+slow!
        for i in xrange(1, N):
            idx = i - 1
            y[i] = y[idx] * mfac + gauss[idx]

        if time_it:
            print(time.time()-t1)

        if array:
            return (y, t)
        else:
            raise NotImplementedError()

    # use slow python implementation for the time being
    # TODO: provide optimized C/weave implementation if possible

    def OU_generator_weave1(self, dt, tau, sigma, y0, t_start=0.0, t_stop=1000.0, time_it=False):
        """
        Generates an Orstein Ulbeck process using the forward euler method. The function returns
        an AnalogSignal object.

        OU_generator_weave1, as opposed to OU_generator, uses scipy.weave
        and is thus much faster.

        Inputs:
            dt      - the time resolution in milliseconds of th signal
            tau     - the correlation time in milliseconds
            sigma   - std dev of the process
            y0      - initial value of the process, at t_start
            t_start - start time in milliseconds
            t_stop  - end time in milliseconds
            array   - if True, the functions returns the tuple (y,t)
                      where y and t are the OU signal and the time bins, respectively,
                      and are both numpy arrays.

        Examples:
            >> stgen.OU_generator_weave1(0.1, 2, 3, 0, 0, 10000)

        See also:
            OU_generator
        """
        import scipy.weave

        import time

        if time_it:
            t1 = time.time()

        t = numpy.arange(t_start, t_stop, dt)
        N = len(t)
        y = numpy.zeros(N, float)
        y[0] = y0
        fac = dt / tau
        gauss = fac * y0 + numpy.sqrt(
            2 * fac) * sigma * self.rng.standard_normal(N - 1)

        # python loop... bad+slow!
        #for i in xrange(1,len(t)):
        # y[i] = y[i-1]+dt/tau*(y0-y[i-1])+numpy.sqrt(2*dt/tau)*sigma*numpy.ran
        # dom.normal()
        # use weave instead
        code = """

        double f = 1.0-fac;

        for(int i=1;i<Ny[0];i++) {
          y(i) = y(i-1)*f + gauss(i-1);
        }
        """

        scipy.weave.inline(code, ['y', 'gauss', 'fac'],
                     type_converters=scipy.weave.converters.blitz)

        if time_it:
            print('Elapsed ', time.time() - t1, ' seconds.')

        if array:
            return (y, t)
        else:
            raise NotImplementedError()

    OU_generator = _OU_generator_python2

    # TODO: optimized inhomogeneous OU generator


# TODO: have a array generator with spatio-temporal correlations

# TODO fix shotnoise stuff below  ... and write tests

# Operations on spike trains

def _gen_g_add(spikes, tau, q, t, eps=1.0e-8):

    #spikes = poisson_generator(rate,t[-1])

    gd_s = numpy.zeros(numpy.shape(t), 'float')

    dt = t[1] - t[0]

    # time of vanishing significance
    vs_t = -tau * log(eps / q)
    kern = q * numpy.exp(-numpy.arange(0.0, vs_t, dt) / tau)

    vs_idx = len(kern)

    idx = numpy.clip(numpy.searchsorted(t, spikes), 0, len(t) - 1)
    idx2 = numpy.clip(idx + vs_idx, 0, len(gd_s))
    idx3 = idx2 - idx

    for i in xrange(len(idx)):

        gd_s[idx[i]:idx2[i]] += kern[0:idx3[i]]

    return gd_s
