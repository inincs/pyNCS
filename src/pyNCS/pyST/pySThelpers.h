
int* extract_cython(int *x, int N, int d, int *a, int *a2, int *r2, int *r, int *y) {
  for (int i=0; i<N; ++i) {
    for (int j=0; j<d; ++j) {
      y[i]+=((x[i] & (a2[j]))>>a[j])*r2[j];
    }
  }

  return y;
}

int* construct_cython(int *x, int N, int d, int *a, int *a2, int *r2, int *r, int *y) {
  for (int i=0; i<N; ++i) {
    for (int j=0; j<d; ++j) {
      y[i]+=((x[i] & (r2[j]))>>r[j])*a2[j];
    }
  }

  return y;
}
