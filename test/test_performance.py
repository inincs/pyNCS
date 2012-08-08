from pyST import *

#stas=STcsMon[0]
#addrPhysicalExtractDecode(stas,[3073]*10000000);

ad=np.random.randint(0, 10000, 10000000)
tm=np.random.ranf(10000000)
t=np.column_stack([ad, tm])
st=SpikeList(t, range(0,20000))


