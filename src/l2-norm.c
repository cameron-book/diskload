#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "diskload.h"
#include <gsl/gsl_integration.h>

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"

#define MICRO "\xc2\xb5"

LoveNumbers* love;
double alpha = 0.1;

double compare_u (double theta, void *params) {
  int cutoff = *(int *) params;
  int output = ((int *) params)[1];
  double w = 1.0;
  double u, v, g;
  double uC, vC, gC;
  
  diskload_truncated( alpha, Uncompensated, theta, w, cutoff, love, &DefaultEarthModel, &u, &v, &g );

  diskload_hypergeometric( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel, &uC, &vC, &gC );

  if (output == 0)
    return (u - uC)*(u-uC);
  if (output == 1)
    return (v - vC)*(v-vC);
  if (output == 2)
    return (g - gC)*(g-gC);    
}

int main( void ) {
  diskload_initialize();
  
  love = diskload_read_love_numbers("../REF_6371_loading_love_numbers_0_40000.txt");
  
  assert( love != NULL );
  assert( love->degrees == 40001 );

  assert( love->h != NULL );  
  assert( love->l != NULL );
  assert( love->k != NULL );  

  diskload_extrapolate_love_numbers(love, 4000000 );

  FILE* table = fopen("table-timings.tex", "w");
  fprintf( table, "\\begin{table}[ht]\n" );
  fprintf( table, "\\centering\n" );  
  fprintf( table, "\\begin{tabular}{lS[table-format=3.2]}\n" );
  fprintf( table, "method & {time (ms)} \\\\\n" );
  fprintf( table, "\\hline\n" );  
 
  double range = 0.2;
  int N = 1000;
  int i;
  double u, v, g;
  double uC, vC, gC;
  double uErr, vErr, gErr;
  DiskLoadError e;

  FILE* f = fopen("l2-norm.csv", "w");
  fprintf( f, "cutoff U V G\n" );

  int cutoff;

  gsl_integration_workspace *integration_workspace = gsl_integration_workspace_alloc (10000);  
  
  //for( cutoff=20000; cutoff<=400000; cutoff += 30000 ) {
  for( cutoff=10000; cutoff<=4000000; cutoff = cutoff*3/2 ) {  
    double lower_limit = 0.50*alpha;	/* lower limit a */
    double upper_limit = 2.00*alpha;	/* upper limit b */
    double abs_error = 1.0e-11;	/* to avoid round-off problems */
    double rel_error = 1.0e-11;	/* the result will usually be much better */
    double result;		/* the result from the integration */
    double error;			/* the estimated error from the integration */
    
    gsl_function My_function;
    int params[2] = { cutoff, 0 };
    void *params_ptr = params;
    My_function.function = &compare_u;
    My_function.params = params_ptr;
    
    double pts[3];
    int npts;
    
    pts[0] = lower_limit;
    pts[1] = upper_limit;
    npts = 2;    

    double u, v, g;
    gsl_integration_qagp (&My_function, pts, npts,
                          abs_error, rel_error, 1000,
                          integration_workspace, &u,
                          &error);

    params[1] = 1;
    gsl_integration_qagp (&My_function, pts, npts,
                          abs_error, rel_error, 1000,
                          integration_workspace, &v,
                          &error);

    params[1] = 2;
    gsl_integration_qagp (&My_function, pts, npts,
                          abs_error, rel_error, 1000,
                          integration_workspace, &g,
                          &error);        

    fprintf( f, "%d %e %e %e\n", cutoff, sqrt(u), sqrt(v), sqrt(g));
    printf( "% 8d %.9f %.9f %.9f\n", cutoff, sqrt(u), sqrt(v), sqrt(g));
  }  
  
  return 0;
}
