# include <R.h>
# include <Rinternals.h>
# include <Rmath.h>
# include <math.h>

void pearson(int* size, double* x, double* y, double* result) {
  double mx, my, xt, yt;
  double xx=0.0, yy=0.0, sxx=0.0, syy=0.0, sxy=0.0; 
  // compute the means
  for(int i=0; i < *size; i++){
    xx += x[i];
    yy += y[i];
  } 
  mx = xx / *size;
  my = yy / *size;
  
  for(int i=0; i < *size; i++){
    xt = x[i] - mx;
    yt = y[i] - my;
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  *result = sxy / (sqrt(sxx*syy));
}
