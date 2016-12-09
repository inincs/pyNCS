from pyNCS.pyST import *
from pyNCS.pyST.STas import load_stas_from_csv, addrLogicalExtract, addrLogicalConstruct, addrPhysicalConstruct, addrPhysicalExtract, addrBuildHashTable
import unittest
import numpy as np
import copy

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        self.stasStim_sac    , self.stasMon_sac     = load_stas_from_csv('chipfiles/sac.csv')
        self.stasStim_ifslwta, self.stasMon_ifslwta = load_stas_from_csv('chipfiles/ifslwta.csv')
        self.stasStim_ifmem  , self.stasMon_ifmem   = load_stas_from_csv('chipfiles/ifmem.csv')
        self.stasSeq         , self.stasMon_linear  = load_stas_from_csv('chipfiles/linear.csv')
        self.stasStim_if2dwta, self.stasMon_if2dwta = load_stas_from_csv('chipfiles/if2dwta.csv')
        self.stasStimRetina  , self.stasRetina      = load_stas_from_csv('chipfiles/tmpdiff64.csv')




        #-------------------------------------stim---- Monitor ----- Unused -- Mon_map ----------------------------
        self.STcsMon=channelAddressing(channelBits=[14,15,16],\
                        stasList=[\
                        self.stasRetina,\
                        self.stasMon_ifslwta,\
                        self.stasMon_ifmem,\
                        self.stasMon_if2dwta,\
                        self.stasSeq,\
                        self.stasStim_ifslwta,\
                        self.stasStim_ifmem,\
                        self.stasStim_if2dwta\
                        ],\
                        )
        #-------------------------------------stim---- Stim_map ---- Unused -- Unused ----------------------------
        self.STcsSeq=channelAddressing(channelBits=[14,15,16],\
                        stasList=[\
                        self.stasRetina,\
                        self.stasStim_ifslwta,\
                        self.stasStim_ifmem,\
                        self.stasStim_if2dwta,\
                        self.stasSeq,\
                        self.stasStim_ifslwta,\
                        self.stasStim_ifmem,\
                        self.stasStim_if2dwta\
                        ],\
                        )
        setDefaultSeqChannelAddress(self.STcsSeq)
        setDefaultMonChannelAddress(self.STcsMon)
        self.tmp_files=[]


        addrHR_out=[range(10),2]
        addr_out=self.STcsMon.addrLogicalConstruct({5:addrHR_out})[5]
        stout=SpikeList(spikes=[],id_list=addr_out)
        time=np.arange(0,50)
        st_out=[]
        for i in addr_out:
            st_out.append((i,STCreate.inh_poisson_generator(t=time,rate=200*np.sin(time*0.01),t_stop=time[-1])))
            stout[i]=st_out[-1][1]

        addrHR=[range(10),5,1]
        addr=self.STcsMon.addrLogicalConstruct({0:addrHR})[0]
        stin=SpikeList(spikes=[],id_list=addr)
        st=[]
        for i in addr:
            st.append((i,STCreate.inh_poisson_generator(t=time,rate=200*np.sin(time*0.01),t_stop=time[-1])))
            stin[i]=st[-1][1]

        self.ch_events={0:stin,1:stout}

    def testEvents_empty(self):
        evs = events()
        evs.add_adtmev([[100,100]])
        evs.empty()

    def testEvents_copy(self):
        evs = events()
        evs.add_adtmev([[100,100],[100,100]])
        events(evs)


    def testStas(self):
        a=addrPhysicalExtract(\
                self.stasStim_ifslwta,addrPhysicalConstruct(\
                self.stasStim_ifslwta,addrLogicalExtract(\
                self.stasStim_ifslwta,addrLogicalConstruct(\
                self.stasStim_ifslwta,[range(5),[15]]\
                ))))
        self.assertTrue(np.all(a==np.array([range(5),[15]*5])))

    def testIsValidAddress(self):
        #Channel 0
        addrHR=[range(10),5,1]
        addrHR_full=[range(10),[5]*10,[1]*10]
        addr=self.STcsMon.addrLogicalConstruct({0:addrHR})[0]        
        self.assertTrue(np.all(addrHR_full==self.STcsMon.addrLogicalExtract({0:addr})[0]))

        #Check that exception is raised when non iterable is given
        self.failUnlessRaises(AssertionError, self.STcsMon.addrLogicalConstruct,{0:1})
        self.failUnlessRaises(AssertionError, self.STcsMon.addrLogicalConstruct,{0:[range(10),5]})

    def testSTCS(self):
        addrHR=[range(10),5,1]
        addr=addrLogicalConstruct(self.STcsMon[0],addrHR)
        addr_out=\
                self.STcsMon.addrLogicalConstruct(\
                self.STcsMon.addrPhysicalExtract(\
                self.STcsMon.addrPhysicalConstruct(\
                self.STcsMon.addrLogicalExtract(\
                {0:addr}\
                ))))[0]
        self.assertTrue(np.all(addr==addr_out))

    def testSTCSChannel3(self):
        #Test 2D chip stim
        addrHR=[range(30,60),2,2]
        addr=addrLogicalConstruct(self.STcsSeq[3],addrHR)
        addr_out=\
                self.STcsSeq.addrLogicalConstruct(\
                self.STcsSeq.addrPhysicalExtract(\
                self.STcsSeq.addrPhysicalConstruct(\
                self.STcsSeq.addrLogicalExtract(\
                {3:addr}\
                ))))[3]
        self.assertTrue(np.all(addr==addr_out))

        #Test 2D chip mon
        addrHR=[range(30,60),2]
        addr=addrLogicalConstruct(self.STcsMon[3],addrHR)
        addr_out=\
                self.STcsMon.addrLogicalConstruct(\
                self.STcsMon.addrPhysicalExtract(\
                self.STcsMon.addrPhysicalConstruct(\
                self.STcsMon.addrLogicalExtract(\
                {3:addr}\
                ))))[3]
        self.assertTrue(np.all(addr==addr_out))
#
#    def test_event_timeslice(self):
#        import pyST
#        import numpy as np
#        a=pyST.STas.events()
#        x=np.arange(0, 500000,5)
#        t=np.arange(100000)*10
#        a.add_adtmev(zip(x,t))
#        sum_t=0
#        for i,j in a.iter_by_timeslice(300):
#            sum_t+=j
#        self.assertTrue(sum_t==t.sum())


#    def testTimeWrapAround(self):
#        addrHR=[range(10),5,1]
#        addr=self.STcsMon.addrLogicalConstruct({0:addrHR})[0]
#        stin=SpikeList(spikes=[],id_list=addr)
#        time=np.arange(0,100)
#        st=[]
#        for i in addr:
#            st.append((i,STCreate.inh_poisson_generator(t=time,rate=100*np.sin(time*0.01),t_stop=100.0)))
#            stin[i]=st[-1][1]
#
#        #Check that SpikeList conserves addresses
#        for element in stin.id_list():
#            self.assert_(element in addr)
#        #Check that Import Export Conserves events
#        for event in st:                
#                events=np.column_stack([event[1].spike_times,np.ones_like(event[1].spike_times)*event[0]])
#                for i in events:
#                    self.assert_( i in stin.raw_data())
#        
#        events=self.STcsSeq.exportAER({0:stin})
#        #Test here: isi, spike conservation
#        aerout=self.STcsSeq.exportAER({0:stin},format='a',isi=False) #returns events        
#        aerout_wrapped=self.STcsSeq.exportAER({0:stin},format='a',isi=False) #returns events        
#        aerout_wrapped.set_tm(aerout.get_tm()%90000)
#        events_out=self.STcsSeq.importAER(aerout_wrapped,format='a',isi=False) #returns channelEvents
#        print aerout.get_tm()
#        for i in events_out.get_tm(0):
#            print i
#            self.assert_(i in aerout.get_tm())


    def testSTCSChannel2(self):
        #Test 2D chip stim
        addrHR=[range(10,20),2,1]
        addr=addrLogicalConstruct(self.STcsSeq[2],addrHR)
        addr_out=\
                self.STcsSeq.addrLogicalConstruct(\
                self.STcsSeq.addrPhysicalExtract(\
                self.STcsSeq.addrPhysicalConstruct(\
                self.STcsSeq.addrLogicalExtract(\
                {2:addr}\
                ))))[2]
        self.assertTrue(np.all(addr==addr_out))

    def testStasSAC(self):
        #Test 2D chip stim
        addrHR=[range(10,20),range(10,15),0]
        addr=addrLogicalConstruct(self.stasStim_sac,addrHR)
        addr_out=\
                self.stasStim_sac.addrLogicalConstruct(\
                self.stasStim_sac.addrPhysicalExtract(\
                self.stasStim_sac.addrPhysicalConstruct(\
                self.stasStim_sac.addrLogicalExtract(\
                addr\
                ))))
        self.assertTrue(np.all(addr==addr_out))

    def testExportImportAER_file(self):
        addrHR=[range(10),5,1]
        addr=self.STcsMon.addrLogicalConstruct({0:addrHR})[0]
        stin=SpikeList(spikes=[],id_list=addr)
        time=np.arange(0,100)
        st=[]
        for i in addr:
            st.append((i,STCreate.inh_poisson_generator(t=time,rate=100*np.sin(time*0.01),t_stop=100.0)))
            stin[i]=st[-1][1]

        #Check that SpikeList conserves addresses
        for element in stin.id_list():
            self.assert_(element in addr)
        #Check that Import Export Conserves events
        for event in st:                
                events=np.column_stack([event[1].spike_times,np.ones_like(event[1].spike_times)*event[0]])
                for i in events:
                    self.assert_( i in stin.raw_data())
        
        events=self.STcsSeq.exportAER({0:stin})
        #Test here: isi, spike conservation
        self.STcsSeq.exportAER({0:stin},filename='aerout',format='a')
        self.tmp_files.append('aerout')
        events_out=self.STcsSeq.importAER(self.tmp_files[0],format='a',isi=True)
        stevents_out=self.STcsSeq.generateST(events_out)

        #Check that Import Export Conserves addresses
        for element in stevents_out[0].select_ids("cell.mean_rate()>0"):
            self.assert_(element in addr)
        
        #Check that Import Export Conserves Events
        for event in st:                
                events=np.column_stack([event[1].spike_times,np.ones_like(event[1].spike_times)*event[0]])
                for i in events:
                    self.assert_( i in stevents_out[0].raw_data())

    def testImportExportEventConservation(self):
        #Create Spike Trains
        events=self.STcsSeq.exportAER(self.ch_events,isi=True)
        events_imported=self.STcsSeq.importAER(events,isi=True)
        #
        self.assert_(self.ch_events[0].raw_data().shape[0]==events_imported[0].get_nev())
        self.assert_(self.ch_events[1].raw_data().shape[0]==events_imported[1].get_nev())

    def testExportOneEvent(self):
        ch_events=SpikeList([(1,817.)],[1])
        events=self.STcsSeq.exportAER(ch_events,isi=True)
        self.assert_(events.get_nev()==1)

    def testEmptyGenerateST(self):
        ch_events=channelEvents(atype='l')
        ch_events.add_adtmch(0,[0],[0])
        self.STcsSeq.generateST(ch_events)
        
    def testNormalizeAER_event_conservation(self):
        events=self.STcsSeq.exportAER(self.ch_events,format='a',isi=True)
        events_imported=self.STcsSeq.importAER(events,isi=True,format='a')

        stnorm=self.STcsSeq.normalizeAER(events_imported)
        self.assert_(self.ch_events[1].raw_data().shape[0]==stnorm[1].get_nev())

        events=self.STcsSeq.exportAER(self.ch_events,format='t',isi=False)
        events_imported=self.STcsSeq.importAER(events,isi=False,format='t')
#        events_imported[1][:,-1]=+100
        #
        stnorm=self.STcsSeq.normalizeAER(events_imported)
        self.assert_(self.ch_events[1].raw_data().shape[0]==stnorm[1].get_nev())

    def testNormalizeAER(self):

        events=self.STcsSeq.exportAER(self.ch_events)
        events_imported=self.STcsSeq.importAER(events,isi=True)

        stnorm=self.STcsSeq.normalizeAER(events_imported)
        for i in stnorm.iterkeys():
            self.assert_(stnorm[i].get_nev()==events_imported[i].get_nev())

        st=self.STcsSeq.generateST(events_imported,normalize=True)
        self.assert_(np.any([st[0].t_start==0,st[1].t_start==0]))

    def testNoEventsExport(self):
        addrHR=[range(10),5,1]
        addr=self.STcsMon.addrLogicalConstruct({0:addrHR})[0]
        stin=SpikeList(spikes=[],id_list=addr)
        time=np.arange(0,100)
        events=self.STcsSeq.exportAER({0:stin})
        self.assertEqual(events.get_nev(),0)

    def testPhysicalLogical(self):
        addrHR=[range(10),5,1]
        addr=self.STcsMon[0].addrPhysicalConstruct(addrHR)
        addrLog=self.STcsMon[0].addrLogicalConstruct(addrHR)
        addrLog_fast=self.STcsMon[0].addrPhysicalLogical(addr)
        for i in addrLog:
            self.assert_( i in addrLog_fast)


    def testPhysicalLogicalSTcsMon(self):    
        addrHR=[range(10),5,1]
        addr=self.STcsMon.addrPhysicalConstruct({0:addrHR})
        addrLog=self.STcsMon.addrLogicalConstruct({0:addrHR})[0]
        addrLog_fast=self.STcsMon.addrPhysicalLogical(addr)[0]
        for i in addrLog:
            self.assert_( i in addrLog_fast)

    def testRawOutputDecode(self):
        #Export an AER stream
        addrHR=[range(10),5,1]
        addr=self.STcsMon.addrLogicalConstruct({0:addrHR})[0]
        stin=SpikeList(spikes=[],id_list=addr)
        time=np.arange(0,100)
        st=[]
        for i in addr:
            st.append((i,STCreate.inh_poisson_generator(t=time,rate=100*np.sin(time*0.01),t_stop=100.0)))
            stin[i]=st[-1][1]


        
        evs = self.STcsSeq.exportAER({0:stin}, isi = False)
        ch_events = self.STcsSeq.extract(evs)
        ch_events = channelEvents(ch_events, atype='p')
        raw_out = self.STcsSeq.rawoutput_from_chevents(ch_events, normalize=False)
        raw_out.decode_all_channels()

        #Check that running this function two times does no harm
        raw_out.decode_all_channels()
        
        tout= np.array(raw_out[0].raw_data())
        tout= np.sort(tout, axis=0).flatten()

        tin = np.array(stin.raw_data())
        tin = np.sort(tin, axis=0).flatten()

        for i,x in enumerate(tin): 
            self.assertAlmostEqual( x, tout[i], delta=1e-3)

        #Check that empty spikelists have complete id_lists
        ## Feature removed in February 2014
        #self.assertEqual(len(raw_out[1].mean_rates()), 128*32)
            
    def testSpikeTrain__time_offset(self):        
        sl = STCreate.poisson_generator(rate=100)
        sl2 = sl.copy()
        sl.time_offset()
        self.assertTrue(np.all(sl.spike_times == sl2.spike_times))
        sl.time_offset(100)
        for i,a in enumerate(sl.spike_times):
            self.assertAlmostEquals(a-100, sl2.spike_times[i], 3)
            
    def testSpikeList__time_offset(self):
        addrHR=[range(10),5,1]
        addr=self.STcsMon.addrLogicalConstruct({0:addrHR})[0]
        sl=SpikeList(spikes=[],id_list=addr)
        for i in addr:            
            sl[i]=STCreate.poisson_generator(rate=10)
        sl2=sl.copy()
        sl.time_offset()
        for i,a in enumerate(sl.raw_data()):
            for c in range(2): 
                self.assertAlmostEquals(a[c], sl2.raw_data()[i][c], 3)
        sl.time_offset(100)
        for i,a in enumerate(sl.raw_data()):
            a[0]-=100
            for c in range(2): 
                self.assertAlmostEquals(a[c], sl2.raw_data()[i][c], 3)


    def testHashTable(self):
        addrHR=[range(64)]
        addrBuildHashTable(self.STcsMon[1])
        addrPhys=self.STcsMon.addrPhysicalConstruct({1:addrHR})
        addrLog=self.STcsMon.addrLogicalConstruct({1:addrHR})[1]
        addrLog_nocheck=self.STcsMon.addrPhysicalLogical(addrPhys)[1]
        addrHR_nocheck=self.STcsMon.addrPhysicalExtract(addrPhys)
        for i in addrLog:
            self.assert_( i in addrLog_nocheck)

    




    def tearDown(self):
        for i in self.tmp_files:
            #os.remove(i)
            pass




if __name__ == '__main__':
    unittest.main()


