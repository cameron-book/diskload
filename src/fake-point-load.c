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

#ifndef MIN_THETA
#define MIN_THETA 10.0
#endif

#ifndef MAX_THETA
#define MAX_THETA 10.0
#endif

#ifndef POSTFIX
#define POSTFIX ""
#endif

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

  diskload_extrapolate_love_numbers(love, 4000000 );

  double range = 0.2;
  int N = 1000;
  int i;
  double w = 1.0;
  double theta;
  double u, v, g;
  DiskLoadError e;

  double alpha = 0.001;
  double min_theta = alpha * MIN_THETA;
  double max_theta = alpha * MAX_THETA;

  double base = exp( log( max_theta / min_theta ) / (N-1) );

  FILE* f = fopen("fake-point-load" POSTFIX ".csv", "w");
  fprintf( f, "alpha theta alpha/theta U V G U40K V40K G40K Upoint Vpoint Gpoint\n" );
  const double Degrees = M_PI / 180.0;
  
  for(i = 0; i<N; i++ ) {
    theta = pow(base, i) * min_theta;

    w = 1.0;
    
    fprintf(f, "%e %e %e ", alpha, theta, theta/alpha);
    e = diskload_hypergeometric( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                            &u, &v, &g );

    fprintf(f, "%e %e %e ", u, v, g);

    printf( "%e ", u );
    double uActual = u;

    w = 1.0;
    
    e = diskload_truncated( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                            &u, &v, &g );

    fprintf(f, "%e %e %e ", u, v, g);

    
    double radius = 2.0 * 1000.0 * DefaultEarthModel.radius * sin(alpha*Degrees / 2.0);
    double area = (M_PI * radius * radius);
    w = area * 4.0 * DefaultEarthModel.gravity / 1e16;

    e = diskload_point( theta, w, 40000, love, &DefaultEarthModel,
                            &u, &v, &g );

    fprintf(f, "%e %e %e ", u, v, g);
    printf( " vs %e (%e)\n", u, uActual/u );    
    fprintf(f, "\n");
  }

  
  return 0;
}
