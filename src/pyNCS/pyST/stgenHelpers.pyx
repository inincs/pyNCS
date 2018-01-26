import cython
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern double* OU_generator_cython(double *y, double *gauss, double fac, int N)

def OU_generator_helper(double[:] y, double[:] gauss, fac, N):
    OU_generator_cython(&y[0], &gauss[0], fac, N);
    return y

