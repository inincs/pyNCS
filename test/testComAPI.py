#!/bin/python
#-----------------------------------------------------------------------------
# File Name : testComAPI.py
# Purpose:
#
# Author: Emre Neftci
#
# Creation Date : 05-06-2013
# Last Modified : Wed 05 Jun 2013 06:47:18 PM PDT
#
# Copyright : (c) 
# Licence : GPLv2
#----------------------------------------------------------------------------- 

from pyNCS.ComAPI import *

class Communicator(BatchCommunicatorBase):
    def run(self, stimulus=None, duration=None, context_manager=None):        
        return stimulus

