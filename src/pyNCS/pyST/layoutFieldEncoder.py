import numpy as np
from scipy import weave

class LayoutFieldEncoder:  # private

    def __init__(self, aspec, nWidth, position=0, pin=''):
        """
        Internal structure used for performing the pin layout permutation 
        Based on the information contained in addrSr 
        Generates the functions "extract" and "construct" to create the final
        bit structure
        """
        #print aspec, nWidth, position, pin
        self.aspec = aspec
        self.aspec2 = 2 ** aspec
        self.nWidth = nWidth
        self.rWidth = np.arange(nWidth)
        self.r2Width = 2 ** self.rWidth
        self.position = position
        self.pin = pin
    
        #Check if transformation is an identity function. If so, just mask and shift
        if np.all(aspec == range(aspec[0],aspec[-1])):
            self.mask = np.sum([2**i for i in aspec]).astype('uint32')
            self.extract = (x&self.mask)>>aspec[0]
            self.construct = (x&self.mask)>>aspec[0]
        else:
            self.extract = lambda x: _extract(x, self.aspec, self.rWidth)
            self.construct = lambda x: _construct(x, self.aspec, self.rWidth)


def _extract(x, a, r):
    '''
    Internal functions for applying the bit permutations defined in the chip file
    '''
    N = int(len(x))
    d = int(r.shape[0])
    r2 = 2 ** r
    a2 = 2 ** a
    y = np.zeros([N], 'int32')

    code = """
           for (int i=0; i<N; ++i) {
               for (int j=0; j<d; ++j) {
                    y[i]+=((x[i] & (a2[j]))>>a[j])*r2[j];
               }
           }
           """
    # compiler keyword only needed on windows with MSVC installed
    ext = weave.inline(code,
                       ['x', 'N', 'd', 'a', 'a2', 'r2', 'r', 'y'],
                       compiler='gcc')
    return y


def _construct(x, a, r):
    '''
    Internal functions for applying the bit permutations defined in the chip file
    '''
    N = int(len(x))
    d = int(r.shape[0])
    r2 = 2 ** r
    a2 = 2 ** a
    y = np.zeros([N], 'int32')

    code = """
           for (int i=0; i<N; ++i) {
               for (int j=0; j<d; ++j) {
                    y[i]+=((x[i] & (r2[j]))>>r[j])*a2[j];
               }
           }
           """
    # compiler keyword only needed on windows with MSVC installed
    ext = weave.inline(code,
                       ['x', 'N', 'd', 'a', 'a2', 'r2', 'r', 'y'],
                       compiler='gcc')
    return y


