#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Sadique Sheik
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
# -*- coding: utf-8 -*-
#Enthougt and Chaco imports
# Work around for Ubuntu 11.10 (are other. distributions affected by the new
# package naming?)

try:
    from enthought.traits.api import HasTraits, Instance, Int, CFloat, Enum, Trait, Callable, Range
    from enthought.traits.ui.api import View, Item, Group, Handler
    from enthought.chaco.api import Plot, ArrayPlotData, jet, Greys
    from enthought.chaco.default_colormaps import *
    from enthought.enable.component_editor import ComponentEditor
    from enthought.traits.ui.menu import Action, CloseAction, MenuBar, Menu
    #Imports for custom marker
    #from enthought.kiva import CompiledPath
except ImportError, e:
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


#def make_custom_marker(size=2):
#    path = CompiledPath()
#    path.move_to(0,size)
#    path.line_to(0,-size)
#    return path
class UpdateEvents:
        def __init__(self, gui, host='localhost', port=50001, channel=0, dims=(128, 0), fps=25):
                self.gui = gui
                self.port = port
                self.host = host
                self.dims = dims
                self.fps = fps
                self.channel = channel
                self.tDuration = 5.
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
                self.neurons = numpy.array([])
                self.times = numpy.array([])

                self.tot_neurons = numpy.array([])
                self.tot_times = numpy.array([])
                #self.z = numpy.zeros(self.dims) #Data being plotted

                self.gui.channel = channel

                #count for updating the rate
                self.clockcount = 0
                self.ratecount = 0

        def fetch(self, *args):
                '''
                The function fetches the data from the server and updates the plot (calles the 'update' function)
                '''
                tDuration = self.gui.tDuration

                cur_neurons = numpy.array([])
                cur_times = numpy.array([])
                # for now i assume that the reading of que is faster than
                # writing.. something more intelligent can go here.
                while True:
                        try:
                                eventsPacket = self.eventsQueue.buffer.get(
                                    block=False)
                                ch_ev = self.stcs.importAER(eventsPacket)
                                if ch_ev.has_key(self.channel):  # Check if there is any data
                                        cur_neurons = append(cur_neurons, ch_ev[self.channel].get_ad(), axis=0,)
                                        cur_times = append(cur_times, ch_ev[self.channel].get_tm())

                        except Queue.Empty, e:
                                break

                self.neurons = cur_neurons
                self.times = cur_times

                old_neurons = self.tot_neurons
                old_times = self.tot_times

                new_times = append(old_times, cur_times * 1e-6)
                new_neurons = append(old_neurons, cur_neurons)

                try:
                        ind = new_times.searchsorted(
                            new_times[-1] - tDuration, side='right')
                        self.tot_neurons = new_neurons[ind:]
                        self.tot_times = new_times[ind:]
                except IndexError, e:
                        pass
                        #print('Warning:Index Error .. (No data)')

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
                        self.gui.plotdata.set_data('neurons', self.tot_neurons)
                        self.gui.plotdata.set_data('times', self.tot_times)
                        self.gui.plot.request_redraw()
                except IndexError, e:
                        print('Warning:Index Error .. (No data)')

                return True

        def stop(self):
                self.eventsQueue.stop()
                print('eventsQueue stopped!')

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
        tDuration = Range(0., 20)
        channel = Enum(range(getDefaultMonChannelAddress().nChannels))

        def _tDuration_default(self):
            return 5.

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
                width=800,
                height=500,
                resizable=True,
                title="Aer1DViewer",
                buttons=['OK'])

        plot_edit_view = View(Group(Item('tDuration'),
                                    #Item('colormap'),
                                    Item('channel')),
                              buttons=['OK', 'Cancel'])

        def __init__(self, dims=(128, 10)):
                super(ImagePlot, self).__init__()
                #z = numpy.zeros(dims)
                self.plotdata = ArrayPlotData(neurons=[0], times=[0])
                plot = Plot(self.plotdata)
                plot.plot(("times", "neurons"),
                          type="scatter",
                          marker="dot",
                          marker_size=1,
                          color='black',)
                self.plot = plot

        def _channel_changed(self):
                print('Switching to channel: %d' % self.channel)
                self.updater.channel = self.channel
                self.updater.tot_neurons = []
                self.updater.tot_times = []
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
