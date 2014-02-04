#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Emre Neftci
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
from __future__ import absolute_import
import os
import numpy as np
import pylab
import matplotlib
from . import pyST_globals
from . import stgen
from .spikes import SpikeList, SpikeTrain, load, merge, merge_spikelists

#Globals
addrIndex = 0
timeIndex = 1
STCreate = stgen.StGen()


def plot_raster(SL, id_list=None, t_start=None, t_stop=None, display=True, id_color=None, kwargs={}):
    ''' Same as spikelist.raster_plot() but with pretty plot options

        id_color is used to customize colors according to neuron address (for example inhibitory vs. excitatory)

        Example:
        >>> id_color = [{'ids':range(124),'color':'blue'},{'ids':range(124,128),'color':'red'}]
        >>> plot_raster(chipout[158],id_color=id_color)
    '''
    #SL._SpikeList__calc_startstop()
    kwargs_default = {'marker': '|', 'markersize': 2, 'color': 'black'}
    kwargs_default.update(kwargs)

    if id_color == None:
        SL.raster_plot(id_list, t_start, t_stop, display, kwargs)
    else:
        h = pylab.axes()
        for i in range(len(id_color)):
            kwargs['color'] = id_color[i]['color']
            SL.id_slice(id_color[i]['ids']).raster_plot(
                t_start=t_start, t_stop=t_stop, display=h, kwargs=kwargs)
        pylab.ylim([SL.id_list().min() - 1, SL.id_list().max() + 1])


STPlotRaster = plot_raster


def composite_plot_superimposed(st_pre, st_post, kwargs={}, colors=['blue', 'k'], id_list=None, t_start=None, t_stop=None, t_start_rate=None, t_stop_rate=None, pre_rate_mult=1.):
    """
    Make a nice Composite plot superimposing two spike trains, a raster plot combined with a vertical rate plot.
    *st_pre* and *st_post* are two SpikeLists. st_post will be plotted above st_pre
    *colors* is a list specifying the respective colors for the bars and spikes for st_pre and st_post
    *kwargs* are the keyword arguments for common plot arguments in the raster plots. Default is {'alpha':0.7,'marker':','} for st_pre and plot_raster defaults for st_post
    *pre_rate_mult* is a float used for plotting the rates of st_pre on a different scale
    other arguments are the same as in composite_plot
    """
    kwargs_pre = dict()
    kwargs_pre.update({'alpha': 0.8, 'marker': ',', 'color': colors[0]})
    kwargs_pre.update(kwargs)
    kwargs_post = dict()
    kwargs_post.update(kwargs)
    kwargs_post.update({'color': colors[1]})

    h1, h2 = composite_plot(
            st_pre,
            kwargs=kwargs_pre,
            kwargs_bar={'alpha': 0.8, 'color': colors[0]},
            id_list=id_list,
            t_start=t_start,
            t_stop=t_stop,
            t_start_rate=t_start_rate,
            t_stop_rate=t_stop_rate,
            )

    twin = False
    if pre_rate_mult != 1.:
        h2t = h2.twiny()
        twin = True
    else:
        h2t = h2

    super_plot = composite_plot(
            st_post,
            kwargs=kwargs_post,
            kwargs_bar={'color': colors[1]},
            display=[h1, h2t],
            id_list=id_list,
            t_start=t_start,
            t_stop=t_stop,
            t_start_rate=t_start_rate,
            t_stop_rate=t_stop_rate)

    if twin:
        pre_lims = h2.get_xlim()
        #h2.set_xlim(pre_lims[0],pre_lims[1]*pre_rate_mult)
        color = colors[0]
        for tl in h2.get_xticklabels():
            tl.set_color(color)
        h2.set_xlabel('Frequency (Hz)', color=color)

    return super_plot


def composite_plot(SL, id_list=None, t_start=None, t_stop=None, t_start_rate=None, t_stop_rate=None, display=True, kwargs={}, kwargs_bar={}):
    """
    Make a nice Composite plot, *i.e.* a raster plot combined with a vertical rate plot.

    The arguments are identical to STPlotRaster, except for display.
    *display* If True is given a new figure is created. If [axS,axR] is given, where axS and axR are pylab.axes objects, then the spike rater and the rate plot are plotted there.
    """
    SL._SpikeList__calc_startstop()

    if id_list == None:
        id_list = SL.id_list()
    if t_start is None:
            t_start = SL.t_start
    if t_stop is None:
            t_stop = SL.t_stop
    if t_start_rate is None:
            t_start_rate = t_start
    if t_stop_rate is None:
            t_stop_rate = t_stop

    if display == True:
        h = pylab.figure()
        axS = pylab.axes([0.12, 0.12, 0.57, 0.8])
        axR = pylab.axes([0.75, 0.12, 0.20, 0.8])
    elif isinstance(display, list):
        axS = display[0]
        axR = display[1]

    STPlotRaster(SL=SL, id_list=id_list, t_start=t_start, t_stop=t_stop,
         display=axS, kwargs=kwargs)

    min_addr = int(np.min(SL.id_list()))
    max_addr = int(np.max(SL.id_list()))
    axS.set_xlim([t_start, t_stop])
    axS.set_yticks([min_addr, (max_addr - min_addr) / 2, max_addr])
    axS.set_xticks(
        [int(t_start), 100 * (int(t_stop - t_start) / 200), int(t_stop)])

    rates = np.array(SL.mean_rates(t_start=t_start_rate, t_stop=t_stop_rate))
    barplot = axR.barh(id_list, rates, linewidth=0, **kwargs_bar)
         #remove space between bars, remove black lines
    pre_lims = axS.set_ylim()
    axR.set_ylim(pre_lims[0], pre_lims[1])

    axR.set_yticks([])
    axR.grid('on')
    #axS.xaxis.grid('on')

    axR.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(2))
    pylab.axes(axR)
    axR.set_xlabel('Frequency (Hz)')

    return axS, axR

STCompositePlot = composite_plot


def composite_plot_movie(SL, time_bin=10, t_start=None, t_stop=None, output="animation.mpg", bounds=(0, 5), fps=10, display=True, maxrate=None, ratebin=None, kwargs={}):
    pylab.ioff()
    from NeuroTools.plotting import get_display
    subplot = get_display(display)
    assert os.system('mencoder') != 32512, "mencoder not found!"
    if t_start is None:
        t_start = SL.t_start
    if t_stop is None:
        t_stop = SL.t_stop
    if maxrate is None:
        maxrate = 100
    if ratebin is None:
        ratebin = 100

    files = []
    im = pylab.figure(**kwargs)
    count = 0
    idx = 0
    axS, axR = STCompositePlot(SL, t_start=0, t_stop=time_bin,
         t_start_rate=0, t_stop_rate=ratebin, kwargs=kwargs)
    axR.hold(False)
    axS.hold(False)
    if t_start != SL.t_start or t_stop != SL.t_stop:
        spk = SL.time_slice(t_start, t_stop)
        raise RuntimeError('Not implemented')

    while (t_start < t_stop):
        STCompositePlot(SL, t_start=0, t_stop=t_start + time_bin, t_start_rate=t_start,
             t_stop_rate=t_start + ratebin, display=[axS, axR], kwargs=kwargs)
        axS.set_xlim([0, t_stop])
        axR.set_xlim([0, maxrate])
        fname = "_tmp_spikes_%05d.png" % count
        pylab.savefig(fname)
        files.append(fname)
        t_start += time_bin
        count += 1
        print('Generated frame {0}'.format(count))
    command = "mencoder 'mf://_tmp_*.png' -mf type=png:fps=%d -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s" % (fps, output)
    os.system(command)
    for fname in files:
        os.remove(fname)


def mapSpikeListAddresses(SL, mapping=None):
    '''
    this function maps the addresses of a spike list into another using the given mapping. Useful for logical to physical translation and vice versa
    SL=original spike list
    mapping=dictionary containing address mapping, If mapping is None, then the addresses will be mapped onto a linear scale, i.e. range(len(SL.id_list()))
    '''
    assert isinstance(SL, SpikeList), "SL must be a NeuroTools SpikeList"
    mapped_SL = SpikeList([], [])
    addr_SL = SL.id_list()
    if mapping == None:
        mapping = dict(zip(np.sort(SL.id_list()), range(len(SL.id_list()))))
    for k, v in mapping.iteritems():
        if k in addr_SL:
            try:
                mapped_SL[v] = SL[k]
            except KeyError:
                pass

    all_others = np.setdiff1d(SL.id_list(), mapping.keys())
    for k in all_others:
        mapped_SL[k] = SL[k]

    return mapped_SL


def ksi(SL, t_stop=800, pow_freq=None, t_bin=10.):
    """
    Kuramoto Synchronization Index (KSI) (Kuaramoto 1984) for measuring phase coherency with respect to the strongest    freqeuncy of the population activity.

    Returns ksi in the range (0,1). 0 meaning no cohrency and 1 meaning entirely coherent

    Inputs:
        SL      - SpikeList of duration of at least t_stop
        t_stop  - the end of the SpikeTrain (in ms)
        t_bin   - time bin for computing population activity and spectrum
        pow_freq- if given, this frequency is taken as the strongest frequency

    """
    #Determine the strongest frequency
    N = len(SL)

    if pow_freq == None:
        psth = SL.time_slice(
            t_start=0, t_stop=t_stop).spike_histogram(t_bin).mean(axis=0)
        t = np.arange(0, t_stop * 1e-3, t_bin * 1e-3)
        sp = abs(np.fft.fft(psth).real)
        freq = abs(np.fft.fftfreq(t.shape[-1], d=t_bin * 1e-3))
        pow_freq = freq[1 + np.argmax(sp.real[1:])]

    l_vs = [[] for i in range(len(SL))]

    for i, k in enumerate(SL.id_list()):
        tm = SL[k].spike_times / 1000  # In seconds
        l_vs[i] = (tm + .5 / pow_freq) % (1. / pow_freq) - .5 / pow_freq

    vs = np.hstack(l_vs) * pow_freq * 2 * np.pi  # Flatten

    return abs(np.sum(np.exp(1j * vs)) / len(vs)), pow_freq
