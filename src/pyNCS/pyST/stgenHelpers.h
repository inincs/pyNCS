
double* OU_generator_cython(double *y, double *gauss, double fac, int N) {
  double f = 1.0-fac;
  
  for(int i=1;i<N;i++) {
    y[i] = y[i-1]*f + gauss[i-1];
  }
  
  return y;
}
