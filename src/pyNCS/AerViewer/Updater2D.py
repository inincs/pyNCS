#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Sadique Sheik
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
try:
    import enthought
    new_version = False
except ImportError:
    new_version = True
# -*- coding: utf-8 -*-
#Enthougt and Chaco imports

if not new_version:
    from enthought.traits.api import HasTraits, Instance, Int, CFloat, Enum, Trait, Callable, Range
    from enthought.traits.ui.api import View, Item, Group, Handler
    from enthought.chaco.api import Plot, ArrayPlotData, jet, Greys
    from enthought.chaco.default_colormaps import *
    from enthought.enable.component_editor import ComponentEditor
    from enthought.traits.ui.menu import Action, CloseAction, MenuBar, Menu
else:
    from traits.api import HasTraits, Instance, Int, CFloat, Enum, Trait, Callable, Range
    from traitsui.api import View, Item, Group, Handler
    from chaco.api import Plot, ArrayPlotData, jet, Greys
    from chaco.default_colormaps import *
    from enable.component_editor import ComponentEditor
    from traitsui.menu import Action, CloseAction, MenuBar, Menu


#Other imports
from numpy import exp, linspace, meshgrid, append
import numpy
import Queue

import wx

#Imports for AER monitoring
import pyAex
from pyNCS.pyST import *
import pyNCS, pyNCS.pyST.STas


class UpdateEvents:
    def __init__(self, gui, host='localhost', port=50001, channel=0, dims=(128, 0), fps=25):
        self.gui = gui
        self.port = port
        self.host = host
        self.dims = dims
        self.fps = fps
        self.channel = channel
        self.gui.updater = self
        self.stcs = getDefaultMonChannelAddress()
        pyNCS.pyST.STas.addrBuildHashTable(self.stcs[channel])

        self.eventsQueue = pyAex.aexclient.AEXMonClient(MonChannelAddress=self.stcs,
                                                   channels=[
                                                       self.channel],
                                                   host=self.host,
                                                   port=self.port,
                                                   autostart=True,
                                                   fps=self.fps)
        self.neurons = numpy.ndarray(shape=(0, 2))
        self.times = numpy.array([])
        self.z = numpy.zeros(self.dims)  # Data being plotted

        self.gui.channel = channel

        #count for updating the rate
        self.clockcount = 0
        self.ratecount = 0

    def fetch(self, *args):
        '''
        The function fetches the data from the server and updates the plot (calles the 'update' function)
        '''

        cur_neurons = numpy.ndarray(shape=(0, 2))
        cur_times = numpy.array([])
        self.datadims = 0  # Dimensionality of the data
        # for now i assume that the reading of que is faster than
        # writing.. something more intelligent can go here.
        while True:
            try:
                eventsPacket = self.eventsQueue.buffer.get(block=False)
                eventsPacket_channel = self.stcs.extract(eventsPacket)
                if self.channel in eventsPacket_channel:  # Check if there is any data
                    pad = eventsPacket_channel.get_ad(self.channel)
                    ad = self.stcs[self.channel].addrPhysicalExtract(pad).transpose()
                    if len(ad[0]) >= 2:  # Confirms that the chip actually has two dimensions.
                        self.datadims = 2
                        cur_neurons = append(cur_neurons,
                                        ad[:, 0:2],
                                        axis=0,)
                    elif len(ad[0]) == 1:  # Confirms that the chip is one dimensional.
                        self.datadims = 1
                        cur_neurons = append(cur_neurons,
                                        [[i, 0] for i in ad[:, 0]],  # crappy hack..
                                        axis=0,)
                    cur_times = append(cur_times,
                               eventsPacket_channel.get_tm(self.channel))
            except Queue.Empty, e:
                    break

        self.neurons = cur_neurons
        self.times = cur_times

        #update clock count
        self.clockcount += 1

        self.update()

    def update(self):
        '''
        update() function updates all the components.
        The class can be inherited and then the update function can be overwritten to make custom filters.
        '''
        self.updatePlot()
        self.updateMeanrate()

    def updateMeanrate(self):
        self.ratecount += len(self.times)
        if self.clockcount % int(self.fps / 2) == 0:
                self.gui.meanrate = self.ratecount * 2  # update the display
                self.ratecount = 0  # reset the counter

    def updatePlot(self):
        '''
        updatePlot() function updates the plot
        '''
        try:
            self.z = self.z * self.gui.decay_factor  # Decay factor
            for c in self.neurons:
                if self.datadims == 2:  # for a 2D chip
                    try:
                        self.z[c[0], c[1]] += 1
                    except IndexError, e:
                        # print('Index error in
                        # updating z',c)
                        pass  # Data out of the dimensions of plot
                elif self.datadims == 1:  # for a 1D chip
                    try:
                        self.z[c[0], 0] += 1
                    except IndexError, e:
                        # print('Index error in
                        # updating z',c)
                        pass  # Data out of the dimensions of plot

            self.gui.plotdata.set_data('imagedata', self.z)
            self.gui.plot.request_redraw()
        except IndexError, e:
            print('Warning:Index Error .. (No data)')
        return True

    def stop(self):
        self.eventsQueue.stop()
        print('Updation ended!')

    def __del__(self):
        self.eventsQueue.stop()


class Controller(Handler):

        view = Instance(HasTraits)

        def init(self, info):
            self.view = info.object

        def edit_plot(self, ui_info):
            self.view.configure_traits(view="plot_edit_view")


class ImagePlot(HasTraits):
        plot = Instance(Plot)
        meanrate = CFloat(0.0)

        updater = Instance(UpdateEvents)

        #Plot properties
        decay_factor = Range(0., 1)
        colormap = Enum(color_map_name_dict.keys())

        channel = Enum(range(getDefaultMonChannelAddress().nChannels))

        #Private properties
        _cmap = Trait(Greys, Callable)

        def _colormap_default(self):
            return 'Greys'

        def _decay_factor_default(self):
            return 0.5

        def _meanrate_default(self):
            return 0.

        traits_view = View(
                Item('plot', editor=ComponentEditor(), show_label=False),
                Item('meanrate', label='MeanRate(Hz)'),
                menubar=MenuBar(Menu(Action(name="Edit Plot",
                                                   action="edit_plot"),
                                            CloseAction,
                                            name="File")),
                handler=Controller,
                width=500,
                height=500,
                resizable=True,
                title="Aer2DViewer",
                buttons=['OK'])

        plot_edit_view = View(Group(Item('decay_factor'),
                                    Item('colormap'),
                                    Item('channel')),
                              buttons=['OK', 'Cancel'])

        def __init__(self, dims=(128, 10)):
                super(ImagePlot, self).__init__()
                z = numpy.zeros(dims)
                self.plotdata = ArrayPlotData(imagedata=z)
                plot = Plot(self.plotdata)
                plot.img_plot("imagedata",
                              xbounds=(0, dims[1]),
                              ybounds=(0, dims[0]),
                              colormap=self._cmap
                              )
                self.plot = plot
                self.flag = True

        def _colormap_changed(self):
                self._cmap = color_map_name_dict[self.colormap]
                if hasattr(self, 'plot'):
                        value_range = self.plot.color_mapper.range
                        self.plot.color_mapper = self._cmap(value_range)
                self.plot.request_redraw()

        def _channel_changed(self):
                print('Switching to channel: %d' % self.channel)
                self.updater.channel = self.channel
                self.updater.z = self.updater.z * 0
                try:
                        self.updater.eventsQueue.stop()
                except:
                        pass
                pyNCS.pyST.STas.addrBuildHashTable(self.updater.stcs[self.channel])


                self.updater.eventsQueue = pyAex.aexclient.AEXMonClient(MonChannelAddress=self.updater.stcs,
                                           channels=[
                                               self.channel],
                                           host=self.updater.host,
                                           port=self.updater.port,
                                           autostart=True,
                                           fps=self.updater.fps)
