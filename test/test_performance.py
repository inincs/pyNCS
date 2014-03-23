from expSetup import *
import time
#stas=STcsMon[0]
#addrPhysicalExtractDecode(stas,[3073]*10000000);
#
#ad=np.random.randint(0, 10000, 10000000)
#tm=np.random.ranf(10000000)
#t=np.column_stack([ad, tm])
#st=SpikeList(t, range(0,20000))
#

nsetup = build_setup(setupfile = 'test_performance.xml')

seq_pop = pyNCS.Population('default', 'Default Population') 
seq_pop.populate_all(nsetup, 'seq', 'excitatory')

fps = 25
t0 = time.time()
for i in range(max(fps,0)):
    stim = seq_pop.soma.spiketrains_poisson(rate=50, duration=1000/fps)
    out = nsetup.run(stim)

t1 = time.time()-t0
print t1

print 'Performance is {0} events per second'.format(float(stim[0].raw_data().__len__())/t1*fps)
