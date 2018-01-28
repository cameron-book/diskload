#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "diskload.h"

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"

#define MICRO "\xc2\xb5"

int main( void ) {
  diskload_initialize();
  
  LoveNumbers* love = diskload_read_love_numbers("../REF_6371_loading_love_numbers_0_40000.txt");
  int tests = 0;
  int successes = 0;
  
  assert( love != NULL );
  assert( love->degrees == 40001 );

  assert( love->h != NULL );  
  assert( love->l != NULL );
  assert( love->k != NULL );  

  diskload_extrapolate_love_numbers(love, 400000 );

  FILE* table = fopen("table-timings.tex", "w");
  fprintf( table, "\\begin{table}[ht]\n" );
  fprintf( table, "\\centering\n" );  
  fprintf( table, "\\begin{tabular}{lS[table-format=3.2]}\n" );
  fprintf( table, "method & {time (ms)} \\\\\n" );
  fprintf( table, "\\hline\n" );  
 
  double range = 0.2;
  int N = 1000;
  int i;
  double w = 1.0;
  double alpha, theta;
  double u, v, g;
  double uC, vC, gC;
  double uErr, vErr, gErr;
  DiskLoadError e;

  FILE* f = fopen("sup-norm.csv", "w");
  fprintf( f, "cutoff U V G\n" );

  int cutoff;
  
  for( cutoff=40000; cutoff<=400000; cutoff += 60000 ) {
    double uMax = 0.0;
    double vMax = 0.0;
    double gMax = 0.0;

    double uTotal = 0.0;
    double vTotal = 0.0;
    double gTotal = 0.0;    
    
    for(i = 0; i<N; i++ ) {
      alpha = 0.3;
      theta = alpha;
      while ( ( fabs(theta - alpha) < 1e-2 ) || (fabs(theta) < 1e-2) )
        //while ( fabs(theta - alpha) < 1e-2 )
        theta = alpha * ((double)rand()/(double)(RAND_MAX/2));
      //}
      e = diskload_truncated( alpha, Uncompensated, theta, w, cutoff, love, &DefaultEarthModel, &u, &v, &g );
      
      e = diskload_hypergeometric( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                                   &uC, &vC, &gC );

      uErr = fabs( uC - u );
      vErr = fabs( vC - v );
      gErr = fabs( gC - g );

      uTotal += ( uC - u )*( uC - u );
      vTotal += ( vC - v )*( vC - v );
      gTotal += ( gC - g )*( gC - g );      
      
      if (uErr > uMax) uMax = uErr;
      if (vErr > vMax) vMax = vErr;
      if (gErr > gMax) gMax = gErr;

      if (i % (N/100) == 0)
        printf( "% 8d %.9f %.9f %.9f\r", cutoff, sqrt(uTotal/i),sqrt(vTotal/i),sqrt(gTotal/i));
    }
    
    fprintf( f, "%d %e %e %e\n", cutoff, sqrt(uTotal),sqrt(vTotal),sqrt(gTotal));    
    printf( "% 8d %.9f %.9f %.9f\n", cutoff, sqrt(uTotal/N),sqrt(vTotal/N),sqrt(gTotal/N));    
  }  
  
  return 0;
}
