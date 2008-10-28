#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <R.h>


void window_metric(double *score, long *nc1_, long *nc2_, double *x1, double *x2, long *w_, long *df_)
{
  /* 
  score - output of the scoring function (nc1 by nc2)
     nc - numbers of cols
      x - vectors with traces
      t - vectors with retention times
      w - penalty for retention times
     df - how far to run the calculation from the diagonal
  */

  int i, j, k, ind, nc1 = *nc1_, nc2=*nc2_, df=*df_, w=*w_, n=(w-1)/2;
  double top=0., bot1=0., bot2=0.;
  
  for (i=0; i< nc1-2*n; i++) {
    for (j=0; j< nc2-2*n; j++) {
   	  if(abs(i-j)>df)
	    continue;
      top = 0.;
      bot1 = 0.;
      bot2 = 0.;
      for (k=0; k<w; k++) {
	    top += x1[i+k]*x2[j+k];
	    bot1 += pow(x1[i+k],2.);
	    bot2 += pow(x2[j+k],2.);
      }
      ind = (i+n)+(j+n)*nc1;
      score[ind] = top / sqrt(bot1*bot2);
    }

  }
}


void cos_ndp_lowmem(double *score, long *nr_, long *nc1_, long *nc2_, double *x1, double *x2, double *t1, double *t2, double *D_, long *df_)
{
  /* 
  score - output of the scoring function (nc1 by nc2)
     nr - number of rows
     nc - numbers of cols
      x - matrices with spectra data
      t - vectors with retention times
      D - penalty for retention times
     df - how far to run the calculation from the diagonal
  */

  int i, j, k, hi, lo, ind, nr = *nr_, nc1 = *nc1_, nc2=*nc2_, df=*df_;
  double dot, D=*D_;
  double *x1sq, *x2sq;
  
  /*
  printf("nr=%d nc1=%d nc2=%d\n",nr,nc1,nc2);
  printf("D=%f df=%d\n",D,df);
  
  double *x1sq = Calloc(*nc1,double);
  double *x2sq = Calloc(*nc2,double);
  */
  
  /* 
  double *x1sq = malloc(nc1 * sizeof(double));
  double *x2sq = malloc(nc2 * sizeof(double));
  */
  x1sq = (double *) R_alloc(nc1, sizeof(double));
  x2sq = (double *) R_alloc(nc2, sizeof(double));

  for (i=0; i< nc1; i++) {
    x1sq[i] = 0.;
    for (k=0; k<nr; k++) {
       x1sq[i] += pow(x1[k+i*nr],2.);
    }
  }
  for (i=0; i< nc2; i++) {
    x2sq[i] = 0.;
    for (k=0; k<nr; k++)
       x2sq[i] += pow(x2[k+i*nr],2.);
  }
  
  /*
  for (i=0; i< nc1; i++)
	printf("x1sq[%d]=%f\n",i,x1sq[i]);
  for (i=0; i< nc2; i++)
	printf("x2sq[%d]=%f\n",i,x2sq[i]);
  for (i=0; i< nc1; i++)
    for (j=0; j< nc2; j++)
	  printf("D[%d,%d]=%f\n",i,j,D[i+j*nr]);
  */
		
  for (i=0; i< nc1; i++) {
    lo = i-df;
    if (lo < 0)
      lo = 0;
    hi = i+df;
    if (hi > nc2)
      hi = nc2;
    
    for (j=lo; j< hi; j++) {
      dot = 0.;
      for (k=0; k<nr; k++) {
        dot += x1[k+i*nr]*x2[k+j*nr];
      }
      ind = i+j*nc1;
      score[ind] = dot / sqrt(x1sq[i]*x2sq[j]);
      score[ind] *= exp(-pow( (t1[i]-t2[j])/D ,2.)/2.);
    }

  } 
   
}
   

void muf(double *v, double *x, long *n_)
{

  int i, j, k, n = *n_, count=0;
  double sum=0.;
  
  /*
  printf("%d\n",n);
  for (i=0; i< n; i++)
	printf("v[%d]=%f\n",i,v[i]);
  for (i=0; i< n*(n+1)/2; i++)
	printf("x[%d]=%f\n",i,x[i]);
  */
  
  for (i=0; i< n; i++)
    for (j=i; j<n; j++) {
	   sum=0.;
       for (k=i; k<=j; k++) {
         sum+=v[k];
	     /* printf("i=%d, j=%d, k=%d, sum=%f\n",i,j,k,sum); */
	   }
	   x[count]=sum/sqrt(j-i+1.);
	   /* printf("(denom=%d) x[%d]=%f\n",(j-i+1),count,x[count]); */
	   count+=1;
    }
   
}
   
   
   
   
void cos_ndp_himem(double *score, long *nr_, long *nc1_, long *nc2_, double *x1, double *x2, double *D, long *df_, double *timedf)
{
  /* 
  score - output of the scoring function (nc1 by nc2)
     nr - number of rows
     nc - numbers of cols
      x - matrices with spectra data
      t - vectors with retention times
      D - penalty for retention times
     df - how far to run the calculation from the diagonal
  */

  int i, j, k, hi, lo, ind, nr = *nr_, nc1 = *nc1_, nc2=*nc2_, df=*df_;
  double dot;
  double *x1sq, *x2sq;
  /*
  double *x1sq = malloc(nc1 * sizeof(double));
  double *x2sq = malloc(nc2 * sizeof(double));
  */
  x1sq = (double *) R_alloc(nc1, sizeof(double));
  x2sq = (double *) R_alloc(nc2, sizeof(double));
 

  for (i=0; i< nc1; i++) {
    x1sq[i] = 0.;
    for (k=0; k<nr; k++) {
       x1sq[i] += pow(x1[k+i*nr],2.);
    }
  }
  for (i=0; i< nc2; i++) {
    x2sq[i] = 0.;
    for (k=0; k<nr; k++)
       x2sq[i] += pow(x2[k+i*nr],2.);
  }
  
  for (i=0; i< nc1; i++) {
    lo = i-df;
    if (lo < 0)
      lo = 0;
    hi = i+df;
    if (hi > nc2)
      hi = nc2;
    
    for (j=lo; j< hi; j++) {
      dot = 0.;
      for (k=0; k<nr; k++) {
        dot += x1[k+i*nr]*x2[k+j*nr];
      }
      ind = i+j*nc1;
      score[ind] = dot / sqrt(x1sq[i]*x2sq[j]);
      score[ind] *= exp(-pow( timedf[i+j*nc1]/D[i+j*nc1] ,2.)/2.);
    }

  } 
   
}


void dp( double *D, double *M, double *gap_, int *phi, long *nr_, long *nc_, int *match, int *nmatch) {
  /* 
      D - output score matrix
      M - input similarity matrix
    gap - gap penalty
     tb - traceback matrix
     nr - number of rows
     nc - number of columns
  */

  int i,j, nr=*nr_, nc=*nc_, tb, count=0, done=0;
  double gap=*gap_, cur_min, a, b, c;

  /*
  for(i=0; i<nr; i++) {
    for(j=0; j<nc; j++) {
	  printf("M[%d,%d]=%f\n",i,j,M[i+j*nr]);
    }
  }
  
  printf("D:\n");
  for(i=0; i<=nr; i++) {
    for(j=0; j<=nc; j++) {
	  printf("%6.5f ",D[i+j*(nr+1)]);
    }
    printf("\n");
  }
  printf("phi:\n");
  for(i=0; i<=nr; i++) {
    for(j=0; j<=nc; j++) {
	  printf("%6d (%6d) ",i+j*(nr+1),phi[i+j*(nr+1)]);
    }
    printf("\n");
  }

  return;
  */
  
  for(i=0; i<nr; i++) {
    for(j=0; j<nc; j++) {
	  a = D[i+j*(nr+1)]+M[i+j*nr];
	  b = D[i+(j+1)*(nr+1)] + gap;
	  c = D[i+1+j*(nr+1)] + gap;
	  /* printf("i=%d j=%d a=%f b=%f c=%f\n",i,j,a,b,c); */
	  if ((a <= b) && (a <= c)) {
	     cur_min=a;
		 tb = 0;
	  }
	  if ((b < a) && (b < c)) {
	     cur_min=b;
		 tb = 1;
	  }
	  if ((c < a) && (c < b)) {
	     cur_min=c;
		 tb = 2;
	  }
	  D[(i+1)+(j+1)*(nr+1)] = cur_min;
	  phi[(i+1)+(j+1)*(nr+1)] = tb;
	  /* printf("D[%d,%d] ind=%d cur_min=%f\n",i+1,j+1,(i+1)+(j+1)*(nr+1),cur_min); */
    }
  }
    
  i=nr;
  j=nc;
  
  while((i > -1) && (j > -1)) {
    /* printf("i=%d j=%d ind=%d\n",i,j,i+j*(nr+1)); */
	switch (phi[i+j*(nr+1)]) {	  
	  case 0:
	    i-=1; j-=1; nmatch[0]++;
		match[count]=i;
		match[count+(nr+nc)]=j;
    	count+=1;
		break;
	  case 1:
		i-=1; break;
	  case 2:
		j-=1; break;
	  case 3:
	    done=1; break;
	}
	if (done==1)
	  break;
  }
  
}

