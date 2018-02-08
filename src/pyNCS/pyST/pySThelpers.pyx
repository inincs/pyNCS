import cython
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern int* extract_cython(int *x, int N, int d, int *a, int *a2, int *r2, int *r, int *y)
cdef extern int* construct_cython(int *x, int N, int d, int *a, int *a2, int *r2, int *r, int *y)

def extract_helper(int[:] x, N, d, int[:] a, int[:] a2, int[:] r2, int[:] r, int[:] y):
    extract_cython(&x[0], N, d, &a[0], &a2[0], &r2[0], &r[0], &y[0]);
    return y

def construct_helper(int[:] x, N, d, int[:] a, int[:] a2, int[:] r2, int[:] r, int[:] y):
    construct_cython(&x[0], N, d, &a[0], &a2[0], &r2[0], &r[0], &y[0]);
    return y;

