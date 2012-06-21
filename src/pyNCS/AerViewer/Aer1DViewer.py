#-----------------------------------------------------------------------------
# Purpose:
#
# Author: Sadique Sheik
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
# -*- coding: utf-8 -*-
import wx


class Aer1DViewer(wx.PySimpleApp):
        '''Class to monitor and display aer activity'''
        def __init__(self, channel=0, port=50001, host='localhost', dims=(-4, 130), fps = 25):
                '''Aer2DViewer class is used to display AER activity in two dimentions
                   channel : channel number.
                   host : Aex server.
                   dims : Range of neurons to be displayed.
                   fps : frames per second.
                   Aer2DViewer(channel=0,
                            host='localhost', dims = (64,32), fps = 40)
                '''
                self.channel = channel
                self.host = host
                self.port = int(port)
                self.dims = dims
                self.fps = fps
                wx.PySimpleApp.__init__(self)

        def OnInit(self, *args, **kw):
                return True

        def setup_timer(self, eventsUpdate):
                # Create a new WX timer
                timerId = wx.NewId()
                self.timer = wx.Timer(self, timerId)

                # Register a callback with the timer event
                self.Bind(wx.EVT_TIMER, eventsUpdate.fetch, id=timerId)

                # Start up the timer!  We have to tell it how many milliseconds
                ## to wait between timer events.
                self.timer.Start(1000 / self.fps, wx.TIMER_CONTINUOUS)
                return

        def show(self):
                '''
                This function pops up the AerViewer
                Alias to MainLoop() function
                '''
                from Updater1D import UpdateEvents, ImagePlot
                self.gui = ImagePlot(dims=self.dims)
                self.eventsUpdate = UpdateEvents(self.gui,
                                               host=self.host,
                                               port=self.port,
                                               channel=self.channel,
                                               dims=self.dims,
                                               fps=self.fps)
                self.setup_timer(self.eventsUpdate)
                self.gui.configure_traits()
                self.eventsUpdate.stop()
                #self.MainLoop()
