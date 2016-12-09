#-----------------------------------------------------------------------------
# Purpose: Test pyNCS functions. Currently using pyAex module, but this should change in the future.
#
# Author: Emre Neftci
#
# Copyright : University of Zurich, Giacomo Indiveri, Emre Neftci, Sadique Sheik, Fabio Stefanini
# Licence : GPLv2
#-----------------------------------------------------------------------------
#import pyAex
import pyNCS
import pyNCS.pyST as pyST
import numpy as np
import unittest, warnings, optparse, os, time, sys



def create_default_population(setup,chipname='seq',N=10,*args,**kwargs):
    '''
    Creates a default population
    '''
    test_pops = pyNCS.Population('default', 'Default Population') 
    test_pops.populate_by_number(setup, chipname, 'excitatory', N, *args, **kwargs)
    return test_pops

def evs_loopback(nsetup, sequencer):
    '''
    This function takes a sequencer object (e.g. returned by
    Group.spiketrains_poisson() ) and creates a RawOutput object out of it. This
    is usefull for generating synthetic data and treating is as neuromorphic
    system ouput.
    '''
    evs = nsetup.mon.exportAER(sequencer, isi=False)
    ch_evs = nsetup.mon.extract(evs)
    return nsetup.mon.rawoutput_from_chevents(ch_evs)


class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        import expSetup
        self.nsetup = expSetup.build_setup()

    def testBuildLinearPopulation(self):
        N=10
        #for transition populations + state populations (for inital and test)
        test_pops=create_default_population(self.nsetup,'seq',N)
        addrs=self.nsetup.mon[test_pops.soma.channel].addrLogicalConstruct([range(N)])
        for i,a in enumerate(test_pops.soma.laddr):
            self.assert_(a in addrs)

    def testMonitors(self):
        N=3
        s=create_default_population(self.nsetup,'seq',N)
        t=create_default_population(self.nsetup, 'ifslwta', N)
        c=pyNCS.PConnection(s,t,'excitatory0')
               
        stmon1=pyNCS.monitors.SpikeMonitor(t.soma)        
        self.nsetup.monitors.import_monitors([stmon1])
        input_stim=s.soma.spiketrains_poisson(rate = np.linspace(10,100,N), duration=500)  
        self.nsetup.prepare()
        self.nsetup.chips['ifslwta'].set_parameter('nsynstdw0',.6)       
        out = self.nsetup.run(input_stim)
        r= stmon1.sl.mean_rates()        
        self.assertTrue(np.all(r>0))

    def testSequencers(self):
        N=3
        s=create_default_population(self.nsetup,'seq',N)
        t=create_default_population(self.nsetup,'ifslwta', N)
        c=pyNCS.PConnection(s,t,'excitatory0')

        sequencer = pyNCS.monitors.Monitors()
        mon_in, = sequencer.create(s)
        mon_in.create_spiketrains('poisson', rate = np.linspace(100,2000,N), duration = 1000)
               
        stmon1=pyNCS.monitors.SpikeMonitor(t.soma)        
        self.nsetup.monitors.import_monitors([stmon1])

        self.nsetup.prepare()
        self.nsetup.chips['ifslwta'].set_parameter('nsynstdw0',.7)       

        out = self.nsetup.run(sequencer)
        r= stmon1.sl.mean_rates()        
        self.assertTrue(np.all(r > 0))

    def testSequencers_nsetup(self):
        N=3
        s=create_default_population(self.nsetup,'seq',N)
        t=create_default_population(self.nsetup,'ifslwta', N)
        c=pyNCS.PConnection(s,t,'excitatory0')

        sequencers = self.nsetup.sequencers
        mon_in, = sequencers.create(s)
        mon_in.create_spiketrains('poisson', rate = np.linspace(100,2000,N), duration = 1000)
               
        stmon1=pyNCS.monitors.SpikeMonitor(t.soma)        
        self.nsetup.monitors.import_monitors([stmon1])

        self.nsetup.prepare()
        self.nsetup.chips['ifslwta'].set_parameter('nsynstdw0',.7)       

        out = self.nsetup.run()
        r= stmon1.sl.mean_rates()        
        self.assertTrue(np.all(r > 0))
                
    def testMonitors_from_SpikeList(self):
        from pyNCS import monitors
        from pylab import close
        N=10
        test_pops=create_default_population(self.nsetup,'seq',N)
        st = test_pops.soma.spiketrains_poisson(10)[0]      
        mon = monitors.create_SpikeMonitor_from_SpikeList(st)        
        monitors.MeanRatePlot(mon)
        close('all')
        
    def testPMapping(self):
        N=30
        p=0.5
        #pyAex.MAPVERS=3
        s=create_default_population(self.nsetup, 'seq', N)
        t=create_default_population(self.nsetup, 'ifslwta', N)
        t2=create_default_population(self.nsetup, 'ifslwta', 124-N, offset=N)
        mon = self.nsetup.monitors.import_monitors_otf(t)[0]
        mon_zero = self.nsetup.monitors.import_monitors_otf(t2)[0]
        m=pyNCS.PMapping('')
        m.connect(s.soma,t.synapses['excitatory0'], fashion='random_all2all', fashion_kwargs={'p':p})
        m.connect(s.soma,t.synapses['excitatory0'], fashion='one2one',fashion_kwargs={'p':p}, expand = True)
        P = int(p*127)
        for i in s.soma.paddr:
            for j in t.synapses['excitatory0'].paddr:
                self.assert_([i, j, P] in m.mapping)
        for n in range(len(s.soma.paddr)):
            self.assert_([s.soma.paddr[n], t.synapses['excitatory0'].paddr[n], P] in m.mapping)
        self.nsetup.mapping.merge(m)
        self.nsetup.prepare()
        self.nsetup.chips['ifslwta'].set_parameter('nsynstdw0',.5)
        input_stim=s.soma.spiketrains_poisson(400)
        out = self.nsetup.run(input_stim)
        self.assertTrue(np.all(350>mon.sl.mean_rates()) and np.all(mon.sl.mean_rates()>100))
        self.assertTrue(np.all(mon_zero.sl.mean_rates()<2))

    def testPConnection(self):
        N=30
        p=0.5
        #pyAex.MAPVERS=3
        s=create_default_population(self.nsetup, 'seq', N)        
        t=create_default_population(self.nsetup, 'ifslwta', N)
        t2=create_default_population(self.nsetup, 'ifslwta', 124-N, offset=N)
        mon = self.nsetup.monitors.import_monitors_otf(t)[0]
        mon_zero = self.nsetup.monitors.import_monitors_otf(t2)[0]
        c=pyNCS.PConnection(s,t,'excitatory0','random_all2all',{'p':p})
        m=c.mapping
        P = int(p*127)
        for i in s.soma.paddr:
            for j in t.synapses['excitatory0'].paddr:
                self.assert_([i, j, P] in m.mapping)
        for n in range(len(s.soma.paddr)):
            self.assert_([s.soma.paddr[n], t.synapses['excitatory0'].paddr[n], P] in m.mapping)

        self.nsetup.prepare()            
        input_stim=s.soma.spiketrains_poisson(400)
        self.nsetup.chips['ifslwta'].set_parameter('nsynstdw0',.5)
        out = self.nsetup.run(input_stim)
        self.assertTrue(np.all(350>mon.sl.mean_rates()) and np.all(mon.sl.mean_rates()>100))
        self.assertTrue(np.all(mon_zero.sl.mean_rates()<2))
        
    def testPMappingLarge(self):
        N=124
        p=0.5
        #pyAex.MAPVERS=3
        s=create_default_population(self.nsetup, 'seq', N)
        t=create_default_population(self.nsetup, 'ifslwta', N)
        m=pyNCS.PMapping('')
        M = np.random.randint(0,2,size=(len(s.soma),2))
        m.mapping.extend(np.random.randint(0,50000,size=(500000,2)))       
        for j in xrange(len(s)):   
            print j      
            m.connect(s.soma[j], t.synapses['inhibitory'][(j*2):((j+1)*2)], fashion = 'by_boolean_matrix', fashion_kwargs={'connection': M[[j],:]})        
        

    def testSeqPopulationFunctions(self):
        N=5
        test_pops1=create_default_population(self.nsetup,'seq',N)
        test_pops2=create_default_population(self.nsetup,'seq',2*N,offset=N)
        pop=pyNCS.Population()
        pop.init(self.nsetup, 'seq', 'excitatory')
        pop.union(test_pops1)
        pop.union(test_pops2)
        testaddr = np.concatenate([test_pops1.soma.paddr,test_pops2.soma.paddr])
        for a in testaddr:
            self.assertTrue(a in pop.soma.paddr)
            
    def testIFSLWTAPopulationFunctions(self):
        N=5
        test_pops1=create_default_population(self.nsetup,'ifslwta',N)
        test_pops2=create_default_population(self.nsetup,'ifslwta',2*N,offset=N)
        pop=pyNCS.Population()
        pop.init(self.nsetup, 'ifslwta', 'excitatory')
        pop.union(test_pops1)
        pop.union(test_pops2)
        testaddr = np.concatenate([test_pops1.soma.paddr,test_pops2.soma.paddr])
        for a in testaddr:
            self.assertTrue(a in pop.soma.paddr)

    def testComAPI_RecordableCommunicatorBase(self):
        import pyNCS.api.ComAPI, os
        rec_com = pyNCS.ComAPI.RecordableCommunicatorBase()
        rec_com.run_rec(np.ones([0,2]))
        self.assertTrue(len(rec_com._rec_fns)==2)
        self.assertTrue(os.path.exists(rec_com._rec_fns[0]))
        self.assertTrue(os.path.exists(rec_com._rec_fns[1]))
        fns = rec_com.get_exp_rec() 
        self.assertTrue(len(rec_com._rec_fns)==0)

    def testExperimentTools(self):
        import pyNCS.experimentTools as et
        test_pops=create_default_population(self.nsetup,'seq', 5)
        input_stim=test_pops.soma.spiketrains_poisson(100)
        self.nsetup.run(input_stim)
        et.mksavedir(pre='/tmp/test_et/')
        et.save_rec_files(self.nsetup)

    def testMonitorsEmptyPlot(self):
        #github inincs/pyNCS issue#3
        test_pops1=create_default_population(self.nsetup,'seq', N=5)
        test_pops2=create_default_population(self.nsetup,'seq', N=5, offset=5)
        stmon1=pyNCS.monitors.SpikeMonitor(test_pops1.soma)
        stmon2=pyNCS.monitors.SpikeMonitor(test_pops2.soma)
        self.nsetup.monitors.import_monitors([stmon1,stmon2])
        mon = test_pops1.soma.spiketrains_poisson(rate = 50)
        #Monitor loopback
        rawoutput = evs_loopback(self.nsetup, mon)
        self.nsetup.monitors.populate_monitors(rawoutput)
        pyNCS.monitors.RasterPlot(self.nsetup.monitors.monitors)
        #import pylab
        #pylab.show()
        
    def testEmptyMonitorReturnsCorrectDimensions(self):
        test_pops1=create_default_population(self.nsetup,'seq', N=5)        
        stmon1=pyNCS.monitors.SpikeMonitor(test_pops1.soma)
        stmon1.populate(pyNCS.pyST.SpikeList())   
        self.assertTrue(len(stmon1.sl.mean_rates())==5)
        




    def tearDown(self):
        del self.nsetup

        
if __name__ == '__main__':
    unittest.main()
    

    #to debug
    #suite.debug()
    #or
    #TestSequenceFunctions('testStimulation').debug()


